library(ggplot2)
library(ggiraph)
library(htmlwidgets)
library(dplyr)
library(jsonlite)

# -----------------------------
# USER SETTINGS
# -----------------------------
gene_name <- "GCG"

chimera_colors <- c(
  "hotpink", "salmon", "cyan", "lime",
  "gold", "violet", "dodgerblue",
  "springgreen", "orange"
)

# -----------------------------
# DATA
# -----------------------------
test <- data.frame(
  AA = LETTERS
) %>%
  mutate(index = row_number())

# -----------------------------
# PLOT
# -----------------------------
p <- ggplot(test) +
  geom_text_interactive(
    aes(x = index, y = 1, label = AA, data_id = index),
    size = 6
  ) +
  theme_void()

wgt <- girafe(
  ggobj = p,
  width_svg = 12,
  height_svg = 4,
  options = list(
    opts_sizing(rescale = FALSE),
    opts_hover(css = "fill:black;"),
    opts_selection(type = "none")
  )
)

# -----------------------------
# CSS
# -----------------------------
wgt$styles <- c(
  wgt$styles,
  "
  text {
    user-select: none;
    pointer-events: all;
  }

  .start-residue {
    fill: orange !important;
    font-weight: bold;
  }

  .selection-rect {
    fill: none;
    stroke-width: 2px;
    rx: 6px;
    ry: 6px;
    pointer-events: none;
  }
  "
)

# -----------------------------
# JAVASCRIPT
# -----------------------------
wgt <- htmlwidgets::onRender(
  wgt,
  paste0("
  const GENE   = '", gene_name, "';
  const COLORS = ", jsonlite::toJSON(chimera_colors), ";

  window.commentLines     = [];
  window.aliasLines       = [];
  window.aliasNames       = [];
  window.selectionRects   = [];
  window.selectedResidues = [];
  window.colorIndex       = 0;

  function setStatus(msg) {
    const el = document.getElementById('selection-status');
    if (el) el.textContent = msg || '';
  }

  function clearStartHighlight() {
    document.querySelectorAll('text.start-residue')
      .forEach(el => el.classList.remove('start-residue'));
  }

  function highlightStartResidue(idx) {
    clearStartHighlight();
    const el = document.querySelector('text[data-id=\"' + idx + '\"]');
    if (el) el.classList.add('start-residue');
  }

  function getTextEl(idx) {
    return document.querySelector('text[data-id=\"' + idx + '\"]');
  }

  function getTextLayer() {
    const svg = document.querySelector('svg');
    if (!svg) return null;
    const g = svg.querySelector('g');
    return g || svg;
  }

  // ---- TIGHT RECT DRAWING ----
  function drawSelectionRect(start, end, color) {
    const t1 = getTextEl(start);
    const t2 = getTextEl(end);
    if (!t1 || !t2) return null;

    const b1 = t1.getBBox();
    const b2 = t2.getBBox();

    const padY  = 2;
    const trimX = 12;

    const x = Math.min(b1.x, b2.x) + trimX / 2;
    const y = Math.min(b1.y, b2.y) - padY;
    const w = Math.abs((b2.x + b2.width) - b1.x) - trimX;
    const h = Math.max(b1.height, b2.height) + padY * 2;

    const rect = document.createElementNS(
      'http://www.w3.org/2000/svg',
      'rect'
    );

    rect.setAttribute('x', x);
    rect.setAttribute('y', y);
    rect.setAttribute('width', w);
    rect.setAttribute('height', h);
    rect.setAttribute('rx', 6);
    rect.setAttribute('ry', 6);
    rect.setAttribute('class', 'selection-rect');
    rect.setAttribute('stroke', color);
    rect.setAttribute('fill', 'none');

    getTextLayer().appendChild(rect);
    return rect;
  }

  function refreshBox() {
    const ta = document.getElementById('saved-result');
    if (!ta) return;

    let out = [];
    out = out.concat(window.commentLines);
    out.push('');
    out = out.concat(window.aliasLines);

    if (window.aliasNames.length) {
      out.push('alias manual_pep ' + window.aliasNames.join('; '));
      out.push('manual_pep');
    }

    ta.value = out.join('\\n');
  }

  // -----------------------------
  // LOAD CXC + RESTORE RECTANGLES
  // -----------------------------
  window.loadCXC = function(file) {
    const reader = new FileReader();
    reader.onload = function(e) {
      window.commentLines   = [];
      window.aliasLines     = [];
      window.aliasNames     = [];
      window.selectionRects.forEach(r => r?.remove());
      window.selectionRects = [];
      window.colorIndex     = 0;

      e.target.result.split(/\\r?\\n/).forEach(line => {
        if (line.startsWith('#')) {
          window.commentLines.push(line);

          const m = line.match(/_(\\d+)-(\\d+)$/);
          if (m) {
            const start = +m[1];
            const end   = +m[2];
            const color = COLORS[window.colorIndex % COLORS.length];
            const rect  = drawSelectionRect(start, end, color);
            window.selectionRects.push(rect);
            window.colorIndex++;
          }
        } else if (
          line.startsWith('alias ') &&
          !line.startsWith('alias manual_pep')
        ) {
          window.aliasLines.push(line);
          const m = line.match(/^alias\\s+(\\S+)/);
          if (m) window.aliasNames.push(m[1]);
        }
      });

      refreshBox();
    };
    reader.readAsText(file);
  };

  // -----------------------------
  // UNDO (REPEATABLE)
  // -----------------------------
  window.undoLast = function() {
    if (!window.commentLines.length) return;

    window.commentLines.pop();
    window.aliasLines.pop();
    window.aliasNames.pop();
    window.colorIndex = Math.max(0, window.colorIndex - 1);

    const r = window.selectionRects.pop();
    if (r) r.remove();

    window.selectedResidues = [];
    clearStartHighlight();
    setStatus('');
    refreshBox();
  };

  // -----------------------------
  // CLICK SELECTION
  // -----------------------------
  document.addEventListener('click', function(e) {
    const t = e.target;
    if (!t || !t.hasAttribute('data-id')) return;

    const idx = +t.getAttribute('data-id');
    if (isNaN(idx)) return;

    window.selectedResidues.push(idx);

    if (window.selectedResidues.length === 1) {
      highlightStartResidue(idx);
      setStatus('Selection started at residue ' + idx);
      return;
    }

    const uniq  = [...new Set(window.selectedResidues)].sort((a,b)=>a-b);
    const start = uniq[0];
    const end   = uniq[uniq.length - 1];

    const alias = 'p' + start + '-' + end;
    if (window.aliasNames.includes(alias)) {
      window.selectedResidues = [];
      clearStartHighlight();
      setStatus('');
      return;
    }

    const color = COLORS[window.colorIndex % COLORS.length];
    window.colorIndex++;

    window.commentLines.push('#' + GENE + '_' + start + '-' + end);

    const deleteLabels =
      window.aliasLines.length === 0 ? 'label delete; ' : '';

    window.aliasLines.push(
      'alias ' + alias + ' ' +
      'select clear; ' +
      deleteLabels +
      'select :' + start + '-' + end + '; ' +
      'color sel ' + color + '; ' +
      'select :' + end + '; ' +
      'label sel text pep_' + start + '-' + end + '; ' +
      'label height 4; ' +
      'select clear'
    );

    window.aliasNames.push(alias);

    const rect = drawSelectionRect(start, end, color);
    window.selectionRects.push(rect);

    window.selectedResidues = [];
    clearStartHighlight();
    setStatus('');
    refreshBox();
  });

  // -----------------------------
  // SAVE
  // -----------------------------
  window.downloadCXC = function() {
    const ta = document.getElementById('saved-result');
    if (!ta?.value) return;

    const blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = GENE + '.cxc';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  // -----------------------------
  // UI
  // -----------------------------
  if (!document.getElementById('annotation-box')) {
    const box = document.createElement('div');
    box.id = 'annotation-box';
    box.style.position = 'fixed';
    box.style.top = '10px';
    box.style.left = '10px';
    box.style.width = '260px';
    box.style.zIndex = 9999;
    box.style.fontFamily = 'monospace';

    box.innerHTML =
      '<input type=\"file\" accept=\".cxc\" onchange=\"loadCXC(this.files[0])\" style=\"width:100%; margin-bottom:4px;\" />' +
      '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
      '<textarea id=\"saved-result\" style=\"width:100%; height:240px;\"></textarea>' +
      '<button onclick=\"undoLast()\" style=\"width:100%; margin-top:4px;\">Undo</button>' +
      '<button onclick=\"downloadCXC()\" style=\"width:100%; margin-top:4px;\">Save</button>';

    document.body.appendChild(box);
  }
  ")
)

# -----------------------------
# SAVE HTML
# -----------------------------
saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)







































library(ggplot2)
library(ggiraph)
library(dplyr)
library(htmlwidgets)
library(jsonlite)

gene_name <- "GCG"

chimera_colors <- c(
  "hotpink", "salmon", "cyan", "lime",
  "gold", "violet", "dodgerblue",
  "springgreen", "orange"
)

test <- data.frame(
  AA = LETTERS
) %>%
  mutate(index = row_number())

p <- ggplot(test) +
  geom_text_interactive(
    aes(x = index, y = 1, label = AA, data_id = index),
    size = 6
  ) +
  theme_void()

wgt <- girafe(
  ggobj = p,
  width_svg = 12,
  height_svg = 4,
  options = list(
    opts_sizing(rescale = FALSE),
    opts_hover(css = "fill:black;"),
    opts_selection(type = "none")
  )
)

wgt$styles <- c(
  wgt$styles,
  "
  text { user-select: text; pointer-events: all; }
  .start-residue {
    fill: orange !important;
    font-weight: bold;
  }
  "
)

wgt <- onRender(
  wgt,
  paste0("
const GENE   = '", gene_name, "';
const COLORS = ", toJSON(chimera_colors), ";

window.commentLines   = [];
window.aliasLines     = [];
window.aliasNames     = [];
window.selectionRects = [];
window.colorIndex     = 0;
window.selectedResidues = [];

const N_RES = document.querySelectorAll('text[data-id]').length;

// ---------- helpers ----------
function setStatus(msg) {
  var el = document.getElementById('selection-status');
  if (el) el.textContent = msg || '';
}

function clearStartHighlight() {
  document.querySelectorAll('text.start-residue')
    .forEach(el => el.classList.remove('start-residue'));
}

function highlightStartResidue(idx) {
  clearStartHighlight();
  document.querySelectorAll('text[data-id=\"' + idx + '\"]')
    .forEach(el => el.classList.add('start-residue'));
}

function refreshBox() {
  var ta = document.getElementById('saved-result');
  if (!ta) return;

  var out = [];
  out = out.concat(window.commentLines);
  out.push('');
  out = out.concat(window.aliasLines);

  if (window.aliasNames.length > 0) {
    out.push('alias manual_pep ' + window.aliasNames.join('; '));
    out.push('manual_pep');
  }

  ta.value = out.join('\\n');
}

// ---------- rectangle track ----------
function redrawSelectionTrack() {
  var svg = document.getElementById('selection-track');
  if (!svg) return;

  svg.innerHTML = '';
  var w = svg.getBoundingClientRect().width;
  var h = svg.getAttribute('height');

  window.selectionRects.forEach(function(sel) {
    var x = (sel.start - 1) / N_RES * w;
    var rw = (sel.end - sel.start + 1) / N_RES * w;

    var r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    r.setAttribute('x', x);
    r.setAttribute('y', 2);
    r.setAttribute('width', Math.max(2, rw));
    r.setAttribute('height', h - 4);
    r.setAttribute('fill', sel.color);
    r.setAttribute('opacity', 0.85);
    svg.appendChild(r);
  });
}

// ---------- load existing cxc ----------
window.loadCXC = function(file) {
  var reader = new FileReader();
  reader.onload = function(e) {
    window.commentLines   = [];
    window.aliasLines     = [];
    window.aliasNames     = [];
    window.selectionRects = [];
    window.colorIndex     = 0;

    var lines = e.target.result.split(/\\r?\\n/);

    lines.forEach(function(line) {
      if (line.startsWith('#')) {
        window.commentLines.push(line);
      } else if (line.startsWith('alias ') && !line.startsWith('alias manual_pep')) {
        window.aliasLines.push(line);
        var m = line.match(/^alias\\s+(p\\d+-\\d+)/);
        if (m) {
          window.aliasNames.push(m[1]);
          var mm = m[1].match(/p(\\d+)-(\\d+)/);
          if (mm) {
            window.selectionRects.push({
              start: +mm[1],
              end:   +mm[2],
              color: COLORS[window.colorIndex % COLORS.length]
            });
            window.colorIndex++;
          }
        }
      }
    });

    refreshBox();
    redrawSelectionTrack();
  };
  reader.readAsText(file);
};

// ---------- undo (repeatable) ----------
window.undoLast = function() {
  if (window.commentLines.length === 0) return;

  window.commentLines.pop();
  window.aliasLines.pop();
  window.aliasNames.pop();
  window.selectionRects.pop();
  window.colorIndex = Math.max(0, window.colorIndex - 1);
  window.selectedResidues = [];

  clearStartHighlight();
  setStatus('');
  refreshBox();
  redrawSelectionTrack();
};

// ---------- click handler ----------
document.addEventListener('click', function(e) {
  var t = e.target;
  if (!t || !t.hasAttribute('data-id')) return;

  var idx = parseInt(t.getAttribute('data-id'), 10);
  if (isNaN(idx)) return;

  window.selectedResidues.push(idx);

  // first click
  if (window.selectedResidues.length === 1) {
    highlightStartResidue(idx);
    setStatus('Selection started at residue ' + idx);
    return;
  }

  var uniq = Array.from(new Set(window.selectedResidues)).sort((a,b) => a-b);
  var start = uniq[0];
  var end   = uniq[uniq.length - 1];

  var alias = 'p' + start + '-' + end;
  if (window.aliasNames.includes(alias)) {
    window.selectedResidues = [];
    clearStartHighlight();
    setStatus('');
    return;
  }

  var color = COLORS[window.colorIndex % COLORS.length];
  window.colorIndex++;

  window.commentLines.push('#' + GENE + '_' + start + '-' + end);

  var deleteLabels =
    window.aliasLines.length === 0 ? 'label delete; ' : '';

  var line =
    'alias ' + alias + ' ' +
    'select clear; ' +
    deleteLabels +
    'select :' + start + '-' + end + '; ' +
    'color sel ' + color + '; ' +
    'select :' + end + '; ' +
    'label sel text pep_' + start + '-' + end + '; ' +
    'label height 4; ' +
    'select clear';

  window.aliasNames.push(alias);
  window.aliasLines.push(line);

  window.selectionRects.push({ start, end, color });

  window.selectedResidues = [];
  clearStartHighlight();
  setStatus('');
  refreshBox();
  redrawSelectionTrack();
});

// ---------- save ----------
window.downloadCXC = function() {
  var ta = document.getElementById('saved-result');
  if (!ta || !ta.value) return;

  var blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = GENE + '.cxc';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
};

// ---------- UI ----------
if (!document.getElementById('annotation-box')) {
  var box = document.createElement('div');
  box.id = 'annotation-box';
  box.style.position = 'fixed';
  box.style.top = '10px';
  box.style.left = '10px';
  box.style.width = '260px';
  box.style.zIndex = 9999;
  box.style.fontFamily = 'monospace';

  box.innerHTML =
    '<input type=\"file\" accept=\".cxc\" onchange=\"loadCXC(this.files[0])\" style=\"width:100%; margin-bottom:4px;\" />' +
    '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
    '<textarea id=\"saved-result\" style=\"width:100%; height:200px; resize:vertical;\"></textarea>' +
    '<svg id=\"selection-track\" width=\"100%\" height=\"24\" style=\"margin-top:6px;border:1px solid #ccc;background:#f9f9f9;\"></svg>' +
    '<button onclick=\"undoLast()\" style=\"width:100%; margin-top:6px;\">Undo last selection</button>' +
    '<button onclick=\"downloadCXC()\" style=\"width:100%; margin-top:4px;\">Save ' + GENE + '.cxc</button>';

  document.body.appendChild(box);
}
")
)

saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)
































library(ggplot2)
library(ggiraph)
library(htmlwidgets)

gene_name <- "GCG"

chimera_colors <- c(
  "hotpink", "salmon", "cyan", "lime",
  "gold", "violet", "dodgerblue",
  "springgreen", "orange"
)

test <- data.frame(
  AA = LETTERS
) %>%
  mutate(index = row_number())

p <- ggplot2::ggplot(test) +
  ggiraph::geom_text_interactive(
    aes(x = index, y = 1, label = AA, data_id = index),
    size = 6
  ) +
  theme_void()

wgt <- ggiraph::girafe(
  ggobj = p,
  width_svg = 12,
  height_svg = 4,
  options = list(
    ggiraph::opts_sizing(rescale = FALSE),
    ggiraph::opts_hover(css = "fill:black;"),
    ggiraph::opts_selection(type = "none")
  )
)

wgt$styles <- c(
  wgt$styles,
  "
  text { user-select: text; pointer-events: all; }
  .start-residue {
    fill: orange !important;
    font-weight: bold;
  }
  "
)

wgt <- htmlwidgets::onRender(
  wgt,
  paste0("
  const GENE = '", gene_name, "';
  const COLORS = ", jsonlite::toJSON(desc_colors), ";

  window.commentLines = [];
  window.aliasLines   = [];
  window.aliasNames   = [];
  window.colorIndex   = 0;
  window.selectedResidues = [];

  function setStatus(msg) {
    var el = document.getElementById('selection-status');
    if (el) el.textContent = msg || '';
  }

  function clearStartHighlight() {
    document.querySelectorAll('text.start-residue')
      .forEach(el => el.classList.remove('start-residue'));
  }

  function highlightStartResidue(idx) {
    clearStartHighlight();
    document.querySelectorAll('text[data-id=\"' + idx + '\"]')
      .forEach(el => el.classList.add('start-residue'));
  }

  function refreshBox() {
    var ta = document.getElementById('saved-result');
    if (!ta) return;

    var out = [];
    out = out.concat(window.commentLines);
    out.push('');
    out = out.concat(window.aliasLines);

    if (window.aliasNames.length > 0) {
      out.push('alias manual_pep ' + window.aliasNames.join('; '));
      out.push('manual_pep');
    }

    ta.value = out.join('\\n');
  }

  // ---------- LOAD EXISTING CXC ----------
  window.loadCXC = function(file) {
    var reader = new FileReader();
    reader.onload = function(e) {
      window.commentLines = [];
      window.aliasLines   = [];
      window.aliasNames   = [];
      window.colorIndex   = 0;

      var lines = e.target.result.split(/\\r?\\n/);

      lines.forEach(function(line) {
        if (line.startsWith('#')) {
          window.commentLines.push(line);
        } else if (line.startsWith('alias ') && !line.startsWith('alias manual_pep')) {
          window.aliasLines.push(line);
          var m = line.match(/^alias\\s+(\\S+)/);
          if (m) {
            window.aliasNames.push(m[1]);
            window.colorIndex++;
          }
        }
      });

      refreshBox();
    };
    reader.readAsText(file);
  };

  // ---------- UNDO ----------
  window.undoLast = function() {
    if (window.commentLines.length === 0) return;

    window.commentLines.pop();
    window.aliasLines.pop();
    window.aliasNames.pop();
    window.colorIndex = Math.max(0, window.colorIndex - 1);
    window.selectedResidues = [];

    clearStartHighlight();
    setStatus('');
    refreshBox();
  };

  // ---------- SELECTION HANDLER ----------
  document.addEventListener('click', function(e) {
    var t = e.target;
    if (!t || !t.hasAttribute('data-id')) return;

    var idx = parseInt(t.getAttribute('data-id'), 10);
    if (isNaN(idx)) return;

    window.selectedResidues.push(idx);

    // FIRST CLICK: highlight
    if (window.selectedResidues.length === 1) {
      highlightStartResidue(idx);
      setStatus('Selection started at residue ' + idx);
      return;
    }

    // SECOND CLICK: commit
    var uniq = Array.from(new Set(window.selectedResidues)).sort((a,b) => a-b);
    var start = uniq[0];
    var end   = uniq[uniq.length - 1];

    var alias = 'p' + start + '-' + end;
    if (window.aliasNames.includes(alias)) {
      window.selectedResidues = [];
      clearStartHighlight();
      setStatus('');
      return;
    }

    var color = COLORS[window.colorIndex % COLORS.length];
    window.colorIndex++;

    window.commentLines.push('#' + GENE + '_' + start + '-' + end);

    var deleteLabels =
      window.aliasLines.length === 0 ? 'label delete; ' : '';

    var line =
      'alias ' + alias + ' ' +
      'select clear; ' +
      deleteLabels +
      'select :' + start + '-' + end + '; ' +
      'color sel ' + color + '; ' +
      'select :' + end + '; ' +
      'label sel text pep_' + start + '-' + end + '; ' +
      'label height 4; ' +
      'select clear';

    window.aliasNames.push(alias);
    window.aliasLines.push(line);

    window.selectedResidues = [];
    clearStartHighlight();
    setStatus('');
    refreshBox();
  });

  // ---------- SAVE ----------
  window.downloadCXC = function() {
    var ta = document.getElementById('saved-result');
    if (!ta || !ta.value) return;

    var blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = GENE + '.cxc';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  // ---------- UI ----------
  if (!document.getElementById('annotation-box')) {
    var box = document.createElement('div');
    box.id = 'annotation-box';
    box.style.position = 'fixed';
    box.style.top = '10px';
    box.style.left = '10px';
    box.style.width = '240px';
    box.style.zIndex = 9999;
    box.style.fontFamily = 'monospace';

    box.innerHTML =
      '<input type=\"file\" accept=\".cxc\" onchange=\"loadCXC(this.files[0])\" style=\"width:100%; margin-bottom:4px;\" />' +
      '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
      '<textarea id=\"saved-result\" style=\"width:100%; height:240px; resize:vertical;\" placeholder=\"Loaded + appended CXC content...\"></textarea>' +
      '<button onclick=\"undoLast()\" style=\"width:100%; margin-top:4px;\">Undo last selection</button>' +
      '<button onclick=\"downloadCXC()\" style=\"width:100%; margin-top:4px;\">Save ' + GENE + '.cxc</button>';

    document.body.appendChild(box);
  }
  ")
)

htmlwidgets::saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)







































library(ggplot2)
library(ggiraph)
library(htmlwidgets)
library(jsonlite)

gene_name <- "GCG"

chimera_colors <- c(
  "hotpink", "salmon", "cyan", "lime",
  "gold", "violet", "dodgerblue",
  "springgreen", "orange"
)

test <- data.frame(
  AA = LETTERS
) %>%
  mutate(index = row_number())

p <- ggplot2::ggplot(test) +
  ggiraph::geom_text_interactive(
    aes(x = index, y = 1, label = AA, data_id = index),
    size = 6
  ) +
  theme_void()

wgt <- ggiraph::girafe(
  ggobj = p,
  width_svg = 12,
  height_svg = 4,
  options = list(ggiraph::opts_sizing(rescale = FALSE),
                 ggiraph::opts_hover(css = "fill:black;"),
                 ggiraph::opts_selection(type = "none"))
)

wgt$styles <- c(
  wgt$styles,
  "
  text { user-select: text; pointer-events: all; }
  .start-residue {
    fill: orange !important;
    font-weight: bold;
  }
  "
)

wgt <- htmlwidgets::onRender(
  wgt,
  "
  const GENE = 'TESTGENE';
  const COLORS = ['hotpink','salmon','skyblue','gold','orchid','limegreen'];

  window.commentLines = [];
  window.aliasLines   = [];
  window.aliasNames   = [];
  window.colorIndex   = 0;
  window.selectedResidues = [];

  function setStatus(msg) {
    var el = document.getElementById('selection-status');
    if (el) el.textContent = msg || '';
  }

  function clearStartHighlight() {
    document.querySelectorAll('text.start-residue')
      .forEach(el => el.classList.remove('start-residue'));
  }

  function highlightStartResidue(idx) {
    clearStartHighlight();
    document.querySelectorAll('text[data-id=\"' + idx + '\"]')
      .forEach(el => el.classList.add('start-residue'));
  }

  function refreshBox() {
    var ta = document.getElementById('saved-result');
    if (!ta) return;

    var out = [];
    out = out.concat(window.commentLines);
    out.push('');
    out = out.concat(window.aliasLines);

    if (window.aliasNames.length > 0) {
      out.push('alias manual_pep ' + window.aliasNames.join('; '));
      out.push('manual_pep');
    }

    ta.value = out.join('\\n');
  }
  // ---------- LOAD EXISTING CXC ----------
  window.loadCXC = function(file) {
    var reader = new FileReader();
    reader.onload = function(e) {
      window.commentLines = [];
      window.aliasLines   = [];
      window.aliasNames   = [];
      window.colorIndex   = 0;

      var lines = e.target.result.split(/\\r?\\n/);

      lines.forEach(function(line) {
        if (line.startsWith('#')) {
          window.commentLines.push(line);
        } else if (line.startsWith('alias ') && !line.startsWith('alias manual_pep')) {
          window.aliasLines.push(line);
          var m = line.match(/^alias\\s+(\\S+)/);
          if (m) {
            window.aliasNames.push(m[1]);
            window.colorIndex++;
          }
        }
      });

      refreshBox();
    };
    reader.readAsText(file);
  };
// ---------- UNDO ----------
  window.undoLast = function() {
    if (window.commentLines.length === 0) return;

    window.commentLines.pop();
    window.aliasLines.pop();
    window.aliasNames.pop();
    window.colorIndex = Math.max(0, window.colorIndex - 1);

    refreshBox();
  };

  // -------- CLICK HANDLER --------
  document.addEventListener('click', function (e) {
    var t = e.target;
    if (!t || !t.hasAttribute('data-id')) return;

    var idx = parseInt(t.getAttribute('data-id'), 10);
    if (isNaN(idx)) return;

    window.selectedResidues.push(idx);

    // ---- FIRST CLICK ----
    if (window.selectedResidues.length === 1) {
      highlightStartResidue(idx);
      setStatus('Selection started at residue ' + idx);
      return;
    }

    // ---- SECOND CLICK ----
    var uniq = Array.from(new Set(window.selectedResidues)).sort((a,b) => a-b);
    var start = uniq[0];
    var end   = uniq[uniq.length - 1];

    var alias = 'p' + start + '-' + end;
    if (window.aliasNames.includes(alias)) {
      window.selectedResidues = [];
      clearStartHighlight();
      setStatus('');
      return;
    }

    var color = COLORS[window.colorIndex % COLORS.length];
    window.colorIndex++;

    window.commentLines.push('#' + GENE + '_' + start + '-' + end);

    var deleteLabels =
      window.aliasLines.length === 0 ? 'label delete; ' : '';

    var line =
      'alias ' + alias + ' ' +
      'select clear; ' +
      deleteLabels +
      'select :' + start + '-' + end + '; ' +
      'color sel ' + color + '; ' +
      'select :' + end + '; ' +
      'label sel text pep_' + start + '-' + end + '; ' +
      'label height 4; ' +
      'select clear';

    window.aliasNames.push(alias);
    window.aliasLines.push(line);

    window.selectedResidues = [];
    clearStartHighlight();
    setStatus('');
    refreshBox();
  });
  // ---------- SAVE ----------
  window.downloadCXC = function () {
    var ta = document.getElementById('saved-result');
    if (!ta || !ta.value) return;

    var blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = GENE + '.cxc';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  // -------- UI --------
  if (!document.getElementById('annotation-box')) {
    var box = document.createElement('div');
    box.id = 'annotation-box';
    box.style.position = 'fixed';
    box.style.top = '10px';
    box.style.left = '10px';
    box.style.width = '220px';
    box.style.zIndex = 9999;
    box.style.fontFamily = 'monospace';

    box.innerHTML =
      '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
      '<textarea id=\"saved-result\" style=\"width:100%;height:240px;\"></textarea>';
      '<textarea id=\"saved-result\" ' +

      'style=\"width:100%; height:240px; resize:vertical;\" ' +
      'placeholder=\"Loaded + appended CXC content...\"></textarea>' +

      '<button onclick=\"undoLast()\" ' +
      'style=\"width:100%; margin-top:4px;\">Undo last selection</button>' +

      '<button onclick=\"downloadCXC()\" ' +
      'style=\"width:100%; margin-top:4px;\">Save ' + GENE + '.cxc</button>';


    document.body.appendChild(box);
  }
  "
)

htmlwidgets::saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)





























library(ggplot2)
library(ggiraph)
library(htmlwidgets)
library(jsonlite)

gene_name <- "HPRDDJSJSS"

chimera_colors <- c(
  "hotpink", "salmon", "cyan", "lime",
  "gold", "violet", "dodgerblue",
  "springgreen", "orange"
)

test <- data.frame(
  AA = paste(rep(LETTERS, 20), collapse = "")
)

p <- ggplot(test) +
  geom_text_interactive(
    aes(x = 1, y = 1, label = AA, data_id = AA),
    size = 6
  ) +
  theme_void()

wgt <- girafe(
  ggobj = p,
  width_svg = 4,
  height_svg = 4,
  options = list(opts_selection(type = "none"), ggiraph::opts_hover(css = "fill:black;"))
)

# Make SVG text selectable
wgt$styles <- c(
  wgt$styles,
  "text { user-select: text; pointer-events: all; }"
)

wgt <- htmlwidgets::onRender(
  wgt,
  paste0("
  const GENE = '", input[["gene"]], "';
  const COLORS = ", jsonlite::toJSON(desc_colors), ";

  window.commentLines = [];
  window.aliasLines   = [];
  window.aliasNames   = [];
  window.colorIndex   = 0;

  function refreshBox() {
    var ta = document.getElementById('saved-result');
    if (!ta) return;

    var out = [];
    out = out.concat(window.commentLines);
    out.push('');
    out = out.concat(window.aliasLines);

    if (window.aliasNames.length > 0) {
      out.push('alias manual_pep ' + window.aliasNames.join('; '));
      out.push('manual_pep');
    }

    ta.value = out.join('\\n');
  }

  // ---------- LOAD EXISTING CXC ----------
  window.loadCXC = function(file) {
    var reader = new FileReader();
    reader.onload = function(e) {
      window.commentLines = [];
      window.aliasLines   = [];
      window.aliasNames   = [];
      window.colorIndex   = 0;

      var lines = e.target.result.split(/\\r?\\n/);

      lines.forEach(function(line) {
        if (line.startsWith('#')) {
          window.commentLines.push(line);
        } else if (line.startsWith('alias ') && !line.startsWith('alias manual_pep')) {
          window.aliasLines.push(line);
          var m = line.match(/^alias\\s+(\\S+)/);
          if (m) {
            window.aliasNames.push(m[1]);
            window.colorIndex++;
          }
        }
      });

      refreshBox();
    };
    reader.readAsText(file);
  };

  // ---------- UNDO ----------
  window.undoLast = function() {
    if (window.commentLines.length === 0) return;

    window.commentLines.pop();
    window.aliasLines.pop();
    window.aliasNames.pop();
    window.colorIndex = Math.max(0, window.colorIndex - 1);

    refreshBox();
  };

  // ---------- SELECTION HANDLER ----------
  document.addEventListener('mouseup', function () {
    var sel = window.getSelection();
    if (!sel || sel.rangeCount === 0) return;

    var range = sel.getRangeAt(0);
    if (!range || range.collapsed) return;

    var node = range.startContainer;
    if (!node || node.nodeType !== Node.TEXT_NODE) return;

    var start = range.startOffset + 1;
    var end   = range.endOffset;
    if (end < start) return;

    var alias = 'p' + start + '-' + end;
    if (window.aliasNames.includes(alias)) return;

    var color = COLORS[window.colorIndex % COLORS.length];
    window.colorIndex++;

    window.commentLines.push('#' + GENE + '_' + start + '-' + end);

    var deleteLabels =
      window.aliasLines.length === 0 ? 'label delete; ' : '';

    var line =
      'alias ' + alias + ' ' +
      'select clear; ' +
      deleteLabels +
      'select :' + start + '-' + end + '; ' +
      'color sel ' + color + '; ' +
      'select :' + end + '; ' +
      'label sel text pep_' + start + '-' + end + '; ' +
      'label height 4; ' +
      'select clear';

    window.aliasNames.push(alias);
    window.aliasLines.push(line);

    refreshBox();
  });

  // ---------- SAVE ----------
  window.downloadCXC = function () {
    var ta = document.getElementById('saved-result');
    if (!ta || !ta.value) return;

    var blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = GENE + '.cxc';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  // ---------- UI ----------
  if (!document.getElementById('annotation-box')) {
    var box = document.createElement('div');
    box.id = 'annotation-box';
    box.style.position = 'fixed';
    box.style.top = '10px';
    box.style.left = '10px';
    box.style.width = '200px';
    box.style.height = '360px';
    box.style.zIndex = 9999;
    box.style.fontFamily = 'monospace';

    box.innerHTML =
      '<input type=\"file\" accept=\".cxc\" ' +
      'onchange=\"loadCXC(this.files[0])\" ' +
      'style=\"width:100%; margin-bottom:4px;\" />' +

      '<textarea id=\"saved-result\" ' +
      'style=\"width:100%; height:240px; resize:vertical;\" ' +
      'placeholder=\"Loaded + appended CXC content...\"></textarea>' +

      '<button onclick=\"undoLast()\" ' +
      'style=\"width:100%; margin-top:4px;\">Undo last selection</button>' +

      '<button onclick=\"downloadCXC()\" ' +
      'style=\"width:100%; margin-top:4px;\">Save ' + GENE + '.cxc</button>';

    document.body.appendChild(box);
  }
  ")
)

htmlwidgets::saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)




































library(ggplot2)
library(ggiraph)
library(htmlwidgets)

gene_name <- "GCG"

# Color palette (ChimeraX-friendly names)
chimera_colors <- c(
  "hotpink", "cyan", "orange", "lime",
  "yellow", "violet", "dodgerblue",
  "springgreen", "gold", "salmon"
)

test <- data.frame(
  AA = paste(c("H", "I", "J", "K", "L", "M", "N"), collapse = "")
)

p <- ggplot(test) +
  geom_text_interactive(
    aes(x = 1, y = 1, label = AA, data_id = AA),
    size = 6
  ) +
  theme_void()

wgt <- girafe(
  ggobj = p,
  width_svg = 4,
  height_svg = 4,
  options = list(opts_selection(type = "none"))
)

# Make text selectable
wgt$styles <- c(
  wgt$styles,
  "text { user-select: text; pointer-events: all; }"
)

wgt <- htmlwidgets::onRender(
  wgt,
  paste0("
  const GENE = '", gene_name, "';
  const COLORS = ", jsonlite::toJSON(chimera_colors), ";
  window.cxcBlocks = window.cxcBlocks || [];
  window.colorIndex = 0;

  function refreshBox() {
    var ta = document.getElementById('saved-result');
    if (ta) ta.value = window.cxcBlocks.join('\\n');
  }

  document.addEventListener('mouseup', function () {
    var sel = window.getSelection();
    if (!sel || sel.rangeCount === 0) return;

    var range = sel.getRangeAt(0);
    if (!range || range.collapsed) return;

    var node = range.startContainer;
    if (!node || node.nodeType !== Node.TEXT_NODE) return;

    var start = range.startOffset + 1;
    var end   = range.endOffset;
    if (end < start) return;

    var color = COLORS[window.colorIndex % COLORS.length];
    window.colorIndex++;

    var header = '#' + GENE + '_' + start + '-' + end;

    var cmd =
      'select clear;\\n' +
      'select :' + start + '-' + end + ';\\n' +
      'color sel ' + color + ';\\n' +
      'select :' + end + ';\\n' +
      'label sel text pep_' + start + '-' + end + ';\\n' +
      'label height 4;\\n' +
      'select clear';

    var block = header + '\\n' + cmd;

    if (!window.cxcBlocks.includes(block)) {
      window.cxcBlocks.push(block);
      refreshBox();
      console.log('Added block:\\n' + block);
    }
  });

  window.downloadCXC = function () {
    if (!window.cxcBlocks.length) return;

    var content = window.cxcBlocks.join('\\n');
    var blob = new Blob([content + '\\n'], { type: 'text/plain' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = GENE + '.cxc';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  if (!document.getElementById('annotation-box')) {
    var box = document.createElement('div');
    box.id = 'annotation-box';
    box.style.position = 'fixed';
    box.style.top = '10px';
    box.style.left = '10px';
    box.style.width = '380px';
    box.style.height = '260px';
    box.style.zIndex = 9999;
    box.style.fontFamily = 'monospace';

    box.innerHTML =
      '<textarea id=\"saved-result\" ' +
      'style=\"width:100%; height:210px; resize:vertical;\" ' +
      'placeholder=\"ChimeraX .cxc commands...\"></textarea>' +
      '<button onclick=\"downloadCXC()\" ' +
      'style=\"width:100%; margin-top:4px;\">Save ' + GENE + '.cxc</button>';

    document.body.appendChild(box);
  }
  ")
)

htmlwidgets::saveWidget(
  wgt,
  file = "~/Desktop/test_annotations.html",
  selfcontained = TRUE
)







library(ggplot2)
library(ggiraph)

p <- ggplot(test) +
  geom_text_interactive(
    aes(
      x = 1,
      y = 1,
      label = AA,
      data_id = AA
    ),
    size = 6
  ) +
  theme_void()

wgt <- girafe(
  ggobj = p,
  options = list(
    opts_selection(type = "none")
  )
)

wgt$styles <- c(
  wgt$styles,
  "
  text {
    user-select: text;
    -webkit-user-select: text;
    pointer-events: all;
  }
  "
)

wgt <- htmlwidgets::onRender(
  wgt,
  "
  // global storage
  window.selectedTexts = window.selectedTexts || [];

  function updateUI() {
    // update visible list
    var panel = document.getElementById(\"saved-result\");
    if (panel) {
      panel.innerHTML =
        '<h3>Selected:</h3><ul>' +
        window.selectedTexts.map(t => '<li>' + t + '</li>').join('') +
        '</ul>';
    }

    // highlight text in SVG
    document.querySelectorAll(\"text\").forEach(t => {
      var txt = t.textContent.trim();
      if (window.selectedTexts.includes(txt)) {
        t.style.fill = \"red\";
        t.style.fontWeight = \"bold\";
      }
    });
  }

  document.addEventListener(\"mouseup\", function () {
    var sel = window.getSelection().toString().trim();
    if (!sel) return;

    sel.split(/\\s+/).forEach(function(s) {
      if (s.length && !window.selectedTexts.includes(s)) {
        window.selectedTexts.push(s);
      }
    });

    updateUI();
    console.log(\"Selected:\", window.selectedTexts);
  });
  "
)

htmlwidgets::saveWidget(wgt, "~/Desktop/test.html", selfcontained = TRUE)

cat(
  '<textarea id="saved-result"
          style="
            position: fixed;
            top: 10px;
            left: 10px;
            width: 300px;
            height: 200px;
            font-family: monospace;
            font-size: 12px;
            z-index: 9999;
          "
          placeholder="Annotations will appear here..."></textarea>

<button onclick="downloadAnnotations()"
        style="position: fixed; top: 220px; left: 10px; z-index: 9999;">
  Save annotations
</button>',
  file = "~/Desktop/test.html",
  append = TRUE
)































test <- tibble(AA = paste0(LETTERS, collapse = ""))

p <- ggplot2::ggplot(data = test) +
  ggiraph::geom_text_interactive(
    aes(
      x = 1,
      y = 1,
      label = AA,
      data_id = AA,
      onclick = "
        var sel = window.getSelection ? window.getSelection().toString() : \"\";
        window.mySavedText = sel;
        console.log(\"Saved text:\", sel);
        var el = document.getElementById(\"saved-result\");
        if (el) {
          el.innerHTML = \"<h3>Saved Text:</h3>\" + sel;
        }
      "
    )
  ) +
  theme_void()

wgt <- ggiraph::girafe(ggobj = p,
                       width_svg = 5,
                       height_svg = 5,
                       options = list(
                         ggiraph::opts_sizing(rescale = FALSE),
                         ggiraph::opts_hover(css = "fill:black;"),
                         ggiraph::opts_selection(type = "none")  # IMPORTANT
                       ))

wgt <- htmlwidgets::onRender(
  wgt,
  "
  // Global storage
  window.selectedTexts = window.selectedTexts || [];

  function updateDisplay() {
    var el = document.getElementById(\"saved-result\");
    if (!el) return;
    el.innerHTML =
      '<h3>Selected:</h3><ul>' +
      window.selectedTexts.map(t => '<li>' + t + '</li>').join('') +
      '</ul>';
  }

  document.addEventListener(\"mouseup\", function () {
    var sel = window.getSelection().toString().trim();
    if (sel.length === 0) return;

    // Optional: split multi-word selections
    sel.split(/\\s+/).forEach(function(s) {
      if (s.length === 0) return;
      if (!window.selectedTexts.includes(s)) {
        window.selectedTexts.push(s);
      }
    });

    updateDisplay();
    console.log(\"Accumulated:\", window.selectedTexts);
  });
  "
)

wgt <- htmlwidgets::onRender(
  wgt,
  "
  window.selectedTexts = [];
  document.addEventListener(\"mouseup\", function () {
    var sel = window.getSelection().toString().trim();
    if (!sel) return;
    window.selectedTexts.push(sel);
    document.getElementById(\"saved-result\").innerHTML =
      '<b>Selected:</b> ' + window.selectedTexts.join(', ');
  });
  "
)

htmlwidgets::saveWidget(widget = wgt,
                        file = "~/Desktop/test.html",
                        selfcontained = TRUE)





p <- ggplot2::ggplot(data = tibble(x = 1, y = 1) %>%
                       mutate(link = "https://stacks.stanford.edu/file/sc075gg6264/test.cxc")) +
  ggiraph::geom_point_interactive(
    aes(
      x = x,
      y = y,
      data_id = x,
      tooltip = link,
      onclick = paste0(
        'window.location="chimerax://open/',
        link,
        '"'
      )
    )
  )



"https://stacks.stanford.edu/file/sc075gg6264/test.cxc"

p <- ggplot2::ggplot(data = tibble(x = 1, y = 1) %>%
                       mutate(link = "https://stacks.stanford.edu/file/sc075gg6264/test.cxc")) +
  ggiraph::geom_point_interactive(aes(x=x,
                                      y=y,
                                      tooltip=paste0("<a href='", link, "'>","Open in ChimeraX",
                                                     "</a>\n"),
                                      onclick=paste0('window.open("', link , '")')))










