import { t } from './i18n.js';

// The survey of starting conditions: an image to look at plus the measured
// numbers behind every pixel (who left, when, how fast neighbours separated).

export class MapPlate {
  constructor(canvas, readout, meta, orbits) {
    this.canvas = canvas;
    this.readout = readout;
    this.meta = meta;
    this.orbits = orbits;
    this.current = null;
    this.hover = null;
    this.onPick = null;
    this.g = canvas.getContext('2d');
    canvas.addEventListener('pointermove', (e) => this.move(e));
    canvas.addEventListener('pointerleave', () => {
      this.hover = null;
      this.draw();
    });
    canvas.addEventListener('click', (e) => {
      const hit = this.locate(e);
      this.onPick?.(hit.p1, hit.p2, hit.cell);
    });
  }

  async load() {
    const [img, bin] = await Promise.all([
      new Promise((resolve, reject) => {
        const im = new Image();
        im.onload = () => resolve(im);
        im.onerror = reject;
        im.src = '/data/map.png';
      }),
      fetch('/data/map.bin').then((r) => r.arrayBuffer()),
    ]);
    this.img = img;
    this.data = new Float32Array(bin);
    this.draw();
  }

  locate(e) {
    const r = this.canvas.getBoundingClientRect();
    const { n, pMax, orderlyBelow } = this.meta;
    const fx = Math.min(Math.max((e.clientX - r.left) / r.width, 0), 0.9999);
    const fy = Math.min(Math.max(1 - (e.clientY - r.top) / r.height, 0), 0.9999);
    const i = Math.floor(fx * n), j = Math.floor(fy * n);
    const o = 3 * (j * n + i);
    const escaper = this.data[o], tEnd = this.data[o + 1], growth = this.data[o + 2];
    return {
      p1: fx * pMax,
      p2: fy * pMax,
      cell: { escaper, tEnd, growth, orderly: escaper < 0 && growth < orderlyBelow },
    };
  }

  move(e) {
    const hit = this.locate(e);
    this.hover = hit;
    const { cell } = hit;
    const fate =
      cell.escaper >= 0
        ? t('map.ejected', { body: cell.escaper + 1, time: cell.tEnd.toFixed(1) })
        : t(cell.orderly ? 'map.regular' : 'map.chaotic', { time: this.meta.tMax });
    this.readout.textContent = t('map.hover', { p1: hit.p1.toFixed(3), p2: hit.p2.toFixed(3), fate });
    this.draw();
  }

  setCurrent(point) {
    this.current = point;
    if (this.img) this.draw();
  }

  toCanvas([p1, p2]) {
    const { pMax } = this.meta;
    return [(p1 / pMax) * this.canvas.width, (1 - p2 / pMax) * this.canvas.height];
  }

  draw() {
    const { g, canvas } = this;
    g.imageSmoothingEnabled = false;
    g.drawImage(this.img, 0, 0, canvas.width, canvas.height);
    g.lineWidth = 1;
    g.strokeStyle = 'rgba(255,255,255,0.7)';
    for (const o of this.orbits) {
      const [x, y] = this.toCanvas(o.params);
      g.beginPath();
      g.arc(x, y, 2.5, 0, Math.PI * 2);
      g.stroke();
    }
    if (this.current) {
      const [x, y] = this.toCanvas(this.current);
      g.strokeStyle = '#fff';
      g.lineWidth = 1.5;
      g.beginPath();
      g.moveTo(x - 12, y);
      g.lineTo(x - 4, y);
      g.moveTo(x + 4, y);
      g.lineTo(x + 12, y);
      g.moveTo(x, y - 12);
      g.lineTo(x, y - 4);
      g.moveTo(x, y + 4);
      g.lineTo(x, y + 12);
      g.stroke();
    }
  }
}
