

chimera_version <- system("ls /Applications | grep ChimeraX", intern = TRUE)
chimera_path <- paste0("/Applications/", chimera_version, "/Contents/MacOS/ChimeraX")



system2(
  chimera_path,
  c("--script",
    cxc_file_path
  )
)







p <- ggplot2::ggplot(data = tibble(x = 0:5, y = 0:5) %>%
                       mutate(link = "https://stacks.stanford.edu/file/sc075gg6264//test_std_chimerax.chimerax.txt")) +
  ggiraph::geom_point_interactive(aes(x = x, y = y, data_id = x, tooltip = link, onclick = paste0('window.open("', link , '")')))



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


wgt <- ggiraph::girafe(ggobj = p,
                       width_svg = 5,
                       height_svg = 5,
                       options = list(
                         ggiraph::opts_sizing(rescale = FALSE)
                       ))

htmlwidgets::saveWidget(widget = wgt,
                        file = "~/Desktop/test.html",
                        selfcontained = TRUE)










