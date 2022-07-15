suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
  library(funr)
})



MyMkdir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x, recursive = T)
  }
}

MyPlotsSave <- function(plot,
                        filename = '__tmp',
                        output = './',
                        width = 10,
                        height = 10,
                        dpi = 500) {
  width = if_else(width > 35, 20, width)
  height = if_else(width > 35, 20, height)
  tryCatch({
    # Give me a plot, save it as ggplot object
    ggsave(as.ggplot(plot, scale = 0.95),
           filename = str_c(filename, '.pdf'),
           path = output,
           width = width,
           height = height,
           limitsize = FALSE
           )
    
    ggsave(
      as.ggplot(plot, scale = 0.95),
      filename = str_c(filename, '.png'),
      path = output,
      width = width,
      height = height,
      dpi = dpi)
  }, 
  error = function(e){
    png(file.path(output, str_c(filename, '.png')), 
        width = width * 18, 
        height = height * 18,
        res = 300,
        units = "mm")
    print(plot)
    dev.off()
    
    pdf(file.path(output, str_c(filename, '.pdf')))
    print(plot)
    dev.off()
  })
}


make_summary <- function(script_path,  status='doing'){
  datetime <- Sys.time()
  script <- basename(script_path)
  pwd <- getwd()
  log <- str_c('datetime:', datetime, '\tpwd:', pwd, '\tscript_id:', script, '\tstatus:', status, '\n')
  log_file <- file.path(dirname(script_path), '../logging_summary.txt')
  cat(log, file = log_file, append = T)
}

