import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';

export const BODY_COLORS = ['#6ec1ff', '#ffc14d', '#ff5c6c'];
const TRAIL_CAP = 40000;

const trailVertex = /* glsl */ `
  attribute float stamp;
  uniform float now;
  uniform float span;
  uniform float braidDrop;
  varying float vFade;
  void main() {
    float age = now - stamp;
    vFade = clamp(1.0 - age / span, 0.0, 1.0);
    vec3 p = position;
    p.z -= age * braidDrop;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(p, 1.0);
  }
`;

const trailFragment = /* glsl */ `
  uniform vec3 tint;
  uniform float gain;
  varying float vFade;
  void main() {
    // vFade can dip a hair below zero where it is extrapolated past a segment's end;
    // pow() of that is undefined and one bad pixel poisons the whole bloom pass
    gl_FragColor = vec4(tint * pow(max(vFade, 0.0), 1.7) * gain, 1.0);
  }
`;

function radialTexture(stops) {
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d');
  const grad = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  for (const [at, alpha] of stops) grad.addColorStop(at, `rgba(255,255,255,${alpha})`);
  g.fillStyle = grad;
  g.fillRect(0, 0, 128, 128);
  return new THREE.CanvasTexture(c);
}

function ringTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 128;
  const g = c.getContext('2d');
  g.strokeStyle = '#fff';
  g.lineWidth = 5;
  g.beginPath();
  g.arc(64, 64, 52, 0, Math.PI * 2);
  g.stroke();
  return new THREE.CanvasTexture(c);
}

class Trail {
  constructor(color) {
    this.positions = new Float32Array(TRAIL_CAP * 3);
    this.stamps = new Float32Array(TRAIL_CAP);
    this.count = 0;
    this.dirtyFrom = 0;
    const geo = new THREE.BufferGeometry();
    this.posAttr = new THREE.BufferAttribute(this.positions, 3).setUsage(THREE.DynamicDrawUsage);
    this.stampAttr = new THREE.BufferAttribute(this.stamps, 1).setUsage(THREE.DynamicDrawUsage);
    geo.setAttribute('position', this.posAttr);
    geo.setAttribute('stamp', this.stampAttr);
    geo.setDrawRange(0, 0);
    this.uniforms = {
      now: { value: 0 },
      span: { value: 10 },
      braidDrop: { value: 0 },
      gain: { value: 1.25 },
      tint: { value: new THREE.Color(color) },
    };
    this.line = new THREE.Line(
      geo,
      new THREE.ShaderMaterial({
        uniforms: this.uniforms,
        vertexShader: trailVertex,
        fragmentShader: trailFragment,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        depthTest: false,
        transparent: true,
      }),
    );
    this.line.frustumCulled = false;
  }

  push(t, x, y, z) {
    if (this.count === TRAIL_CAP) {
      const keep = TRAIL_CAP >> 1;
      this.positions.copyWithin(0, (TRAIL_CAP - keep) * 3);
      this.stamps.copyWithin(0, TRAIL_CAP - keep);
      this.count = keep;
      this.dirtyFrom = 0;
    }
    const i = this.count++;
    this.positions[3 * i] = x;
    this.positions[3 * i + 1] = y;
    this.positions[3 * i + 2] = z;
    this.stamps[i] = t;
  }

  flush() {
    if (this.dirtyFrom < this.count) {
      for (const [attr, size] of [[this.posAttr, 3], [this.stampAttr, 1]]) {
        attr.clearUpdateRanges();
        attr.addUpdateRange(this.dirtyFrom * size, (this.count - this.dirtyFrom) * size);
        attr.needsUpdate = true;
      }
      this.dirtyFrom = this.count;
    }
    this.line.geometry.setDrawRange(0, this.count);
  }

  clear() {
    this.count = 0;
    this.dirtyFrom = 0;
    this.line.geometry.setDrawRange(0, 0);
  }
}

export class Stage {
  constructor(container) {
    this.container = container;
    const renderer = (this.renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' }));
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color('#10122e');
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.01, 4000);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(0.6, -3.2, 2.4);

    this.controls = new OrbitControls(this.camera, renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.06;
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 0.35;

    this.composer = new EffectComposer(renderer);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloom = new UnrealBloomPass(new THREE.Vector2(1, 1), 0.6, 0.55, 0.2);
    this.composer.addPass(this.bloom);
    this.composer.addPass(new OutputPass());

    this.scale = 1;
    this.braid = 0;
    this.braidTarget = 0;
    this.span = 10;
    this.inset = 0;

    this.buildStars();
    this.buildDial();
    this.buildBodies();

    new ResizeObserver(() => this.resize()).observe(container);
    this.resize();
  }

  buildStars() {
    const n = 1800;
    const pos = new Float32Array(n * 3);
    const col = new Float32Array(n * 3);
    const c = new THREE.Color();
    for (let i = 0; i < n; i++) {
      const u = Math.random() * 2 - 1, phi = Math.random() * Math.PI * 2;
      const r = 900 + Math.random() * 900, q = Math.sqrt(1 - u * u);
      pos.set([r * q * Math.cos(phi), r * q * Math.sin(phi), r * u], 3 * i);
      const warm = Math.random();
      c.setHSL(warm < 0.7 ? 0.62 : 0.09, 0.2, 0.12 + 0.5 * Math.random() ** 4);
      col.set([c.r, c.g, c.b], 3 * i);
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(col, 3));
    this.stars = new THREE.Points(
      geo,
      new THREE.PointsMaterial({ size: 1.4, sizeAttenuation: false, vertexColors: true, depthWrite: false }),
    );
    this.scene.add(this.stars);
  }

  // A faint engraved dial in the plane the bodies start in, for a sense of scale.
  buildDial() {
    const pts = [];
    const ring = (r, seg = 180) => {
      for (let i = 0; i < seg; i++) {
        const a = (i / seg) * Math.PI * 2, b = ((i + 1) / seg) * Math.PI * 2;
        pts.push(r * Math.cos(a), r * Math.sin(a), 0, r * Math.cos(b), r * Math.sin(b), 0);
      }
    };
    [0.5, 1, 1.5, 2].forEach((r) => ring(r));
    for (let i = 0; i < 72; i++) {
      const a = (i / 72) * Math.PI * 2, r0 = 2, r1 = i % 6 === 0 ? 2.12 : 2.05;
      pts.push(r0 * Math.cos(a), r0 * Math.sin(a), 0, r1 * Math.cos(a), r1 * Math.sin(a), 0);
    }
    pts.push(-2, 0, 0, 2, 0, 0, 0, -2, 0, 0, 2, 0);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
    this.dial = new THREE.LineSegments(
      geo,
      new THREE.LineBasicMaterial({ color: '#ece6d6', transparent: true, opacity: 0.07, depthWrite: false }),
    );
    this.scene.add(this.dial);
  }

  buildBodies() {
    const glowTex = radialTexture([[0, 0.9], [0.1, 0.3], [0.35, 0.05], [1, 0]]);
    const ringTex = ringTexture();
    const sphere = new THREE.SphereGeometry(1, 32, 16);
    this.bodies = [];
    this.trails = [];
    this.paths = [];
    this.ghosts = [];
    BODY_COLORS.forEach((hex) => {
      const core = new THREE.Mesh(sphere, new THREE.MeshBasicMaterial({ color: new THREE.Color(hex).lerp(new THREE.Color('#fff'), 0.3) }));
      const glow = new THREE.Sprite(
        new THREE.SpriteMaterial({ map: glowTex, color: hex, blending: THREE.AdditiveBlending, depthWrite: false, depthTest: false }),
      );
      core.add(glow);
      this.scene.add(core);
      this.bodies.push(core);

      const trail = new Trail(hex);
      this.scene.add(trail.line);
      this.trails.push(trail);

      const path = new THREE.Line(
        new THREE.BufferGeometry(),
        new THREE.LineBasicMaterial({ color: hex, transparent: true, opacity: 0.2, blending: THREE.AdditiveBlending, depthWrite: false }),
      );
      path.frustumCulled = false;
      path.visible = false;
      this.scene.add(path);
      this.paths.push(path);

      const ghost = new THREE.Sprite(
        new THREE.SpriteMaterial({ map: ringTex, color: hex, blending: THREE.AdditiveBlending, depthWrite: false, depthTest: false }),
      );
      ghost.visible = false;
      this.scene.add(ghost);
      this.ghosts.push(ghost);
    });
  }

  // Size everything to a system whose bodies stay within `extent` of the centre.
  frame(extent, masses) {
    this.scale = extent;
    this.dial.scale.setScalar(extent / 1.5);
    this.bodies.forEach((b, i) => {
      const r = 0.022 * extent * Math.cbrt(masses[i]);
      b.scale.setScalar(r);
      b.children[0].scale.setScalar(13);
      this.ghosts[i].scale.setScalar(0.16 * extent);
    });
    const dir = this.camera.position.clone().sub(this.controls.target).normalize();
    this.controls.target.set(0, 0, 0);
    this.camera.position.copy(dir.multiplyScalar(extent * 3.4));
    this.controls.minDistance = extent * 0.3;
    this.controls.maxDistance = extent * 60;
  }

  setPath(samples) {
    this.paths.forEach((path, b) => {
      if (!samples) {
        path.visible = false;
        return;
      }
      const pts = new Float32Array(samples.length * 3);
      samples.forEach((s, i) => pts.set([s[3 * b], s[3 * b + 1], s[3 * b + 2]], 3 * i));
      path.geometry.dispose();
      path.geometry = new THREE.BufferGeometry();
      path.geometry.setAttribute('position', new THREE.BufferAttribute(pts, 3));
      path.visible = true;
    });
  }

  setBodies(s) {
    this.bodies.forEach((b, i) => b.position.set(s[3 * i], s[3 * i + 1], s[3 * i + 2]));
  }

  setGhosts(points) {
    this.ghosts.forEach((g, i) => {
      g.visible = !!points;
      if (points) g.position.set(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
    });
  }

  pushTrail(t, s) {
    this.trails.forEach((tr, i) => tr.push(t, s[3 * i], s[3 * i + 1], s[3 * i + 2]));
  }

  clearTrails() {
    this.trails.forEach((tr) => tr.clear());
  }

  resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    this.renderer.setSize(w, h);
    this.composer.setSize(w, h);
    this.camera.aspect = w / h;
    // keep the system centred in the part of the window the drawer leaves free
    if (this.inset) this.camera.setViewOffset(w, h, this.inset / 2, 0, w, h);
    else this.camera.clearViewOffset();
    this.camera.updateProjectionMatrix();
  }

  setInset(px) {
    this.inset = px;
    this.resize();
  }

  render(now) {
    this.braid += (this.braidTarget - this.braid) * 0.08;
    const drop = (this.braid * 2.6 * this.scale) / this.span;
    this.trails.forEach((tr) => {
      tr.flush();
      tr.uniforms.now.value = now;
      tr.uniforms.span.value = this.span;
      tr.uniforms.braidDrop.value = drop;
    });
    this.dial.material.opacity = 0.07 * (1 - 0.6 * this.braid);
    this.paths.forEach((p) => (p.material.opacity = 0.2 * (1 - this.braid)));
    this.controls.target.z = -1.1 * this.braid * this.scale;
    this.controls.update();
    this.composer.render();
  }
}
