


html_slide_show <- function(svg_directory = "~/huBEC/ehbp_normalBECwFetal/animate2",
                            output_file = "~/Desktop/test.html",
                            title = NULL,
                            frames = sub(".svg", "", list.files(svg_directory)),
                            buttonNames = frames,
                            categories = NULL,
                            columns = "auto") {

  if(is.null(categories)) categories <- rep("none", length(frames))

  doc <- xml2::read_html(system.file("templates/slide_show_template.html", package = "ligandFinder"))

  if(!is.null(title)) {
    titleNode <- xml2::xml_find_first(doc, "//title")
    xml2::xml_text(titleNode) <- title
  }

  scriptNode <- xml2::xml_find_first(doc, "//script")

  buttonNode <- xml2::xml_find_first(doc, "//div[@class='controls2']")

  if(columns > 1) {
    cols <- cut(1:length(frames), columns, labels = FALSE)
  } else {
    cols <- rep(1, length(frames))
  }

  for(x in unique(cols)) {
    xml2::xml_add_child(.x = buttonNode,
                        xml2::read_xml(paste0("<div class='column' id='col", x, "'></div>")),
                        .where = "after")
  }

  uni_cats <- c(unique(categories), "banana")
  cat_index <- 1
  for(f in seq_along(frames)) {
    buttonNode <- xml2::xml_find_first(doc, paste0("//div[@class='controls2']//div[@id='col", cols[f], "']"))

    if(categories[f] == uni_cats[cat_index]) {
      if(!uni_cats[cat_index] == "none") {
        xml2::xml_add_child(.x = buttonNode,
                            xml2::read_xml(paste0("<h2>", uni_cats[cat_index], "</h2>")))
      }
      xml2::xml_add_child(.x = buttonNode,
                          xml2::read_xml(paste0("<button>", buttonNames[f], "</button>")))
      cat_index <- cat_index + 1
    } else {
      xml2::xml_add_child(.x = buttonNode,
                          xml2::read_xml(paste0("<button>", buttonNames[f], "</button>")))
    }
  }



  for(f in seq_along(frames)) {
    svgDoc <- xml2::read_xml(paste0(svg_directory, "/", frames[f], ".svg"))
    xml2::xml_add_child(.x = xml2::xml_find_first(doc, "//div[@class='svg-container']"), svgDoc)
  }

  xml2::write_html(doc, output_file)

}


