
make_html_loader <- function(
    html_dir = "htmls",
    output   = "index.html"
) {
  files <- list.files(
    html_dir,
    pattern = "\\.html$",
    ignore.case = TRUE,
    full.names = TRUE
  )

  stopifnot(length(files) > 0)

  out_dir <- dirname(output)

  # paths RELATIVE to index.html location
  rel_paths <- vapply(
    files,
    function(f) {
      p <- normalizePath(f, winslash = "/", mustWork = TRUE)
      o <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
      sub(paste0("^", o, "/?"), "", p)
    },
    character(1)
  )

  options <- paste0(
    '<option value="', rel_paths, '">',
    sub("\\.html$", "", basename(files), ignore.case = TRUE),
    '</option>',
    collapse = "\n"
  )

  html <- paste0(
    '<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>HTML Loader</title>

<style>
body { margin: 0; font-family: sans-serif; }
#wrap { display: flex; height: 100vh; }
#left { width: 250px; padding: 10px; border-right: 1px solid #ccc; }
#right { flex: 1; overflow: auto; }
iframe { border: none; width: 100%; height: 100%; }
select { font-family: monospace; }
</style>
</head>

<body>
<div id="wrap">
  <div id="left">
    <select id="files" size="25" style="width:100%">
', options, '
    </select>
  </div>
  <div id="right">
    <iframe id="frame"></iframe>
  </div>
</div>

<script>
const sel = document.getElementById("files");
const frame = document.getElementById("frame");

sel.addEventListener("change", () => {
  frame.src = sel.value;
});

if (sel.options.length > 0) {
  sel.selectedIndex = 0;
  frame.src = sel.value;
}
</script>
</body>
</html>'
  )

  writeLines(html, output)
  invisible(output)
}


setwd("~/Desktop")
make_html_loader(
  html_dir = "tmp2",
  output = "tmp2/index.html"
)



