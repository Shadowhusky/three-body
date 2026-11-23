const STRINGS = {
  en: {
    title: 'Three bodies',
    'tab.orbits': 'Orbits',
    'tab.map': 'Map',
    'tab.findings': 'Findings',
    switchTo: '中文',
    'orbits.lede': 'Periodic orbits found by our own search. Each one closes on itself to within rounding error.',
    'orbits.stableOnly': 'Stable only',
    'orbits.count': '{shown} of {total}',
    'orbits.none': 'No orbits loaded',
    'orbits.period': 'period {period}',
    stable: 'stable',
    unstable: 'unstable',
    'group.flat': 'Planar orbits',
    'group.flat.note': 'Zero angular momentum. All three bodies move in one plane.',
    'group.family': 'One family',
    'group.family.note':
      'A continuous family of 3D orbits, each periodic in a rotating frame. Every entry is the same pattern rotating at a different rate; any fraction can be computed on request.',
    'group.spatial': '3D orbits',
    'group.spatial.note': 'Systems with angular momentum whose paths leave the plane.',
    'group.orderly': 'Quasi-periodic orbits',
    'group.orderly.note': 'They never repeat, yet stay regular. A formula fitted to a short span predicts them far ahead.',
    'map.lede':
      'Each point is one choice of initial velocities. Colour shows which body is eventually ejected; ivory marks systems that stay bound and regular. Click anywhere to watch that system.',
    'map.axisX': 'outer bodies’ velocity along the line',
    'map.axisY': 'outer bodies’ velocity across the line',
    'map.hint': 'Move over the map to inspect a system.',
    'map.hover': 'Velocity ({p1}, {p2}): {fate}.',
    'map.ejected': 'body {body} is ejected at time {time}',
    'map.regular': 'bound and regular through time {time}',
    'map.chaotic': 'still bound at time {time}, but chaotic',
    'read.time': 'Time',
    'read.energy': 'Energy error',
    'read.sep': 'Closest separation',
    'read.formula': 'Formula error',
    'read.tape': 'Who passed between the other two',
    pause: 'Pause',
    play: 'Play',
    restart: 'Restart',
    speed: 'Speed',
    trail: 'Trail',
    braid: 'Braid through time',
    showFormula: 'Show formula',
    autoRotate: 'Auto-rotate',
    'note.period': 'Period {period}.',
    'note.stable': 'Linearly stable: small disturbances stay small.',
    'note.unstable': 'Unstable: a small disturbance grows {growth}-fold per period.',
    'note.family': 'Family member: the whole pattern rotates by {turn} of a half-turn every half cycle.',
    'note.quasi':
      'Quasi-periodic: it never repeats. A {terms}-term formula fitted to the first {fitSpan} time units stays within {error} of the true motion out to time {testedTo}. Turn on "Show formula" to compare.',
    'pick.name': 'Initial velocity ({p1}, {p2})',
    'pick.ejected': 'In the survey, body {body} is ejected around time {time}.',
    'pick.regular': 'The survey classifies this system as regular: bound and predictable.',
    'pick.chaotic': 'The survey classifies this system as chaotic: still bound when the survey ended.',
    fate: 'Body {body} is ejected at time {time}. The other two remain bound as a binary.',
    noData: 'No research data loaded yet.',
  },
  zh: {
    title: '三体',
    'tab.orbits': '轨道',
    'tab.map': '图谱',
    'tab.findings': '发现',
    switchTo: 'EN',
    'orbits.lede': '由我们自己的搜索程序找到的周期轨道。每一条都在舍入误差之内首尾闭合。',
    'orbits.stableOnly': '只看稳定的',
    'orbits.count': '{shown} / {total}',
    'orbits.none': '未加载轨道',
    'orbits.period': '周期 {period}',
    stable: '稳定',
    unstable: '不稳定',
    'group.flat': '平面轨道',
    'group.flat.note': '角动量为零，三个天体始终在同一平面内运动。',
    'group.family': '同一个轨道族',
    'group.family.note':
      '一族连续的三维轨道，每一条在旋转参考系中都是周期的。这里的每一项都是同一种图样以不同的速度旋转；任意分数都可以按需算出来。',
    'group.spatial': '三维轨道',
    'group.spatial.note': '带有角动量、轨迹离开平面的系统。',
    'group.orderly': '拟周期轨道',
    'group.orderly.note': '永不重复，却始终规则。用一小段运动拟合出的公式，可以预测很久以后。',
    'map.lede':
      '每个点代表一组初始速度。颜色表示最终哪个天体被抛出；象牙色是保持束缚、运动规则的系统。点击任意位置，观看那个系统。',
    'map.axisX': '外侧两天体沿连线方向的速度',
    'map.axisY': '外侧两天体垂直于连线的速度',
    'map.hint': '把鼠标移到图上，查看某个系统。',
    'map.hover': '速度 ({p1}, {p2})：{fate}。',
    'map.ejected': '天体 {body} 在时间 {time} 被抛出',
    'map.regular': '到时间 {time} 仍然束缚且规则',
    'map.chaotic': '到时间 {time} 仍然束缚，但已混沌',
    'read.time': '时间',
    'read.energy': '能量误差',
    'read.sep': '最近距离',
    'read.formula': '公式误差',
    'read.tape': '谁从另外两个中间穿过',
    pause: '暂停',
    play: '播放',
    restart: '重新开始',
    speed: '速度',
    trail: '轨迹',
    braid: '沿时间编成辫子',
    showFormula: '显示公式',
    autoRotate: '自动旋转',
    'note.period': '周期 {period}。',
    'note.stable': '线性稳定：小扰动不会放大。',
    'note.unstable': '不稳定：小扰动每个周期放大 {growth} 倍。',
    'note.family': '轨道族成员：整个图样每半个周期转过半圈的 {turn}。',
    'note.quasi':
      '拟周期：永不重复。用前 {fitSpan} 个时间单位拟合出的 {terms} 项公式，直到时间 {testedTo} 都与真实运动相差不超过 {error}。打开“显示公式”对比。',
    'pick.name': '初始速度 ({p1}, {p2})',
    'pick.ejected': '在巡查中，天体 {body} 大约在时间 {time} 被抛出。',
    'pick.regular': '巡查把这个系统归为规则：保持束缚，可以预测。',
    'pick.chaotic': '巡查把这个系统归为混沌：巡查结束时仍然束缚。',
    fate: '天体 {body} 在时间 {time} 被抛出，另外两个仍然束缚成双星。',
    noData: '尚未加载研究数据。',
  },
};

const ORBIT_NAMES_ZH = [
  [/^Figure-eight$/, '8 字形'],
  [/^Butterfly (.+)$/, '蝴蝶 $1'],
  [/^Moth (.+)$/, '飞蛾 $1'],
  [/^Yin-yang (.+)$/, '阴阳 $1'],
  [/^Dragonfly$/, '蜻蜓'],
  [/^Goggles$/, '护目镜'],
  [/^Yarn$/, '线团'],
  [/^Bumblebee$/, '熊蜂'],
  [/^Planar orbit (\d+)$/, '平面轨道 $1'],
  [/^3D orbit (\d+)$/, '三维轨道 $1'],
  [/^Turn (.+)$/, '转速 $1'],
  [/^Quasi-periodic orbit (\d+)$/, '拟周期轨道 $1'],
];

const TAGS_ZH = [
  [/^known orbit$/, '已知轨道'],
  [/^found by search$/, '搜索发现'],
  [/^computed from the family$/, '由轨道族算出'],
  [/^not in the 2025 catalogue$/, '未收录于 2025 年目录'],
  [/^catalogue orbit (O\d+)$/, '目录轨道 $1'],
  [/^existence proved$/, '存在性已证明'],
  [/^predicted by a (\d+)-term formula$/, '可用 $1 项公式预测'],
];

function stored() {
  try {
    return localStorage.getItem('lang');
  } catch {
    return null;
  }
}

let lang = new URLSearchParams(location.search).get('lang') || stored() || (navigator.language.startsWith('zh') ? 'zh' : 'en');
if (!STRINGS[lang]) lang = 'en';

export const getLang = () => lang;

export function setLang(next) {
  lang = next;
  try {
    localStorage.setItem('lang', next);
  } catch {
    // private browsing: the choice just won't persist
  }
}

export function t(key, vars = {}) {
  return (STRINGS[lang][key] ?? STRINGS.en[key]).replace(/\{(\w+)\}/g, (_, name) => vars[name]);
}

function translate(text, table) {
  for (const [pattern, replacement] of table) {
    if (pattern.test(text)) return text.replace(pattern, replacement);
  }
  return text;
}

// Orbit names and tags are stored in English in the data; these give the display form.
export function orbitName(name) {
  return lang === 'zh' ? translate(name, ORBIT_NAMES_ZH) : name;
}

export function orbitTag(tag) {
  if (lang !== 'zh' || !tag) return tag;
  return tag
    .split(', ')
    .map((part) => translate(part, TAGS_ZH))
    .join('，');
}

// Fill every element marked data-i18n with the current language.
export function applyStaticText() {
  document.documentElement.lang = lang === 'zh' ? 'zh-CN' : 'en';
  document.title = t('title');
  document.querySelectorAll('[data-i18n]').forEach((el) => (el.textContent = t(el.dataset.i18n)));
}
