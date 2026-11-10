// Exercise the embedded report without browser packages.
// Usage: node descriptor_exploration_ui_test.js /path/to/report.html
'use strict';
const fs = require('fs'), vm = require('vm'), assert = require('assert');
const page = fs.readFileSync(process.argv[2], 'utf8');
const data = JSON.parse(page.match(/<script id="report-data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const script = page.match(/<script>([\s\S]*?)<\/script>/)[1];
data.page_size = 1;
class Element {
  constructor(tag) {
    this.tag = tag;
    this.children = [];
    this.events = {};
    this.value = '';
  }
  append(...nodes) { this.children.push(...nodes); }
  replaceChildren(...nodes) { this.children = nodes; }
  addEventListener(name, callback) { this.events[name] = callback; }
  setAttribute() {}
  scrollIntoView() {}
  getContext() { return new Proxy({}, {get: () => () => {}}); }
  getBoundingClientRect() { return {left: 0, top: 0, width: 650, height: 285}; }
  click() {}
}
const elements = {};
for (const id of ['report-data', 'search', 'evidence', 'sort', 'rows', 'count', 'prev', 'next', 'export', 'inspect']) {
  elements[id] = new Element('div');
}
elements['report-data'].textContent = JSON.stringify(data);
elements.sort.value = 'score';
let captured;
const context = vm.createContext({
  document: {getElementById: id => elements[id], createElement: tag => new Element(tag)},
  Blob,
  URL: {createObjectURL: blob => {captured = blob; return 'blob:test';}, revokeObjectURL() {}},
  setTimeout: callback => callback(),
  console,
});
vm.runInContext(script, context);
const names = () => Array.from(vm.runInContext('selected().map(row => row.name)', context));
const expected = data.features.slice().sort((a, b) => b.score-a.score || a.name.localeCompare(b.name));
assert.deepEqual(names(), expected.map(row => row.name));
assert.equal(elements.rows.children.length, 2);
assert.equal(elements.rows.children[0].children[0].textContent, expected[0].name);
if (expected.length > 1) {
  elements.next.onclick();
  assert.equal(elements.rows.children[0].children[0].textContent, expected[1].name);
  elements.prev.onclick();
}
function change(id, value) {elements[id].value = value; elements[id].events.input();}
change('search', 'no matching descriptor exists');
assert.equal(elements.rows.children.length, 0);
change('search', '');
change('evidence', 'weak');
assert.deepEqual(names(), expected.filter(row => row.evidence === 'weak').map(row => row.name));
change('evidence', '');
change('sort', 'pearson');
const magnitudes = Array.from(vm.runInContext('selected().map(row => Math.abs(row.pearson))', context));
assert(magnitudes.every((value, i) => i === 0 || magnitudes[i-1] >= value));
const detail = elements.rows.children[1].children[0].children[0];
detail.open = true;
detail.events.toggle(); // Exercise curve and interval drawing.
vm.runInContext("inspect(data.features[0].buckets[0].groups[0].members, 'Selected bucket')", context);
const members = data.features[0].buckets[0].groups[0].members;
assert.equal(elements.inspect.children[0].textContent, 'Selected bucket');
assert.equal(elements.inspect.children[3].children.length, Math.min(3, members.length));
elements.inspect.children[2].onclick();
const bucket = captured;
elements.export.onclick();
const ranking = captured;
Promise.all([bucket.text(), ranking.text()]).then(([a, b]) => {
  assert.equal(a.trim().split('\n').length, members.length + 1);
  for (const i of members) assert(a.includes('"'+data.observations.identifiers[i]+'"'));
  assert.equal(b.trim().split('\n').length, data.features.length + 1);
  console.log('Report filtering, sorting, pagination, plotting, examples and complete CSV exports passed');
}).catch(error => {console.error(error); process.exitCode = 1;});
