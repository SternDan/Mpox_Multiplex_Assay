calculateCutofflog10 <- function(data, antigene){
  detCut <- data %>% 
    #    filter(population == "Kupferzell") %>%
    filter(!is.na({{antigene}})) %>%
    dplyr::select({{antigene}})
  detCut <- log10(detCut[[1]])
  detCut_out <- em(detCut,"normal","normal")
  confint <- confint(detCut_out,nb=100,level=.95)
  cut_off <- cutoff(detCut_out)
  return(list("input" = detCut, "out" = detCut_out, "confint" = confint, "cut_off" = cut_off))
}

plotCutoff <- function(input, outname, type_fun, type, width, height){
  type_fun(paste(outname, type, sep = "."), width = width, height = height)
  hist(input$input, 100, F, xlab = "Bioplex (MFI)", ylab = "Density", main = NULL, col = "grey")
  lines(input$out,lwd=1.5,col="red")
  polygon(c(input$cut_off[-1],rev(input$cut_off[-1])),c(0,0,1,1),col=rgb(0,0,1,.2),border=NA)
  abline(v=input$cut_off[-1],lty=2,col="blue")
  abline(v=input$cut_off[1],col="blue")
  dev.off()
}

# Function and plotting for ELISAs: keine log-transformation der Daten
calculateCutoff <- function(data, antigene){
  detCut <- data %>% 
    #    filter(population == "Kupferzell") %>%
    filter(!is.na({{antigene}})) %>%
    dplyr::select({{antigene}})
  detCut <- (detCut[complete.cases(detCut),])
  detCut_out <- em(detCut,"normal","normal")
  confint <- confint(detCut_out,nb=100,level=.95)
  cut_off <- cutoff(detCut_out)
  return(list("input" = detCut, "out" = detCut_out, "cut_off" = cut_off))
}

plotCutoffELISA <- function(input, outname, type_fun, type, width, height, cutSupplier){
  type_fun(paste(outname, type, sep = "."), width = width, height = height)
  hist(input$input, 100, F, xlab = "ELISA", ylab = "Density", main = NULL, col = "grey")
  lines(input$out,lwd=1.5,col="red")
  polygon(c(input$cut_off[-1],rev(input$cut_off[-1])),c(0,0,1,1),col=rgb(0,0,1,.2),border=NA)
  abline(v=input$cut_off[-1],lty=2,col="blue")
  abline(v=input$cut_off[1],col="blue")
  abline(v=cutSupplier,col="green")
  dev.off()
}

plotCutoffEuroimmun <- function(input, outname, type_fun, type, width, height, lower, upper){
  type_fun(paste(outname, type, sep = "."), width = width, height = height)
  hist(input$input, 100, F, xlab = "ELISA", ylab = "Density", main = NULL, col = "grey")
  lines(input$out,lwd=1.5,col="red")
  polygon(c(input$cut_off[-1],rev(input$cut_off[-1])),c(0,0,1,1),col=rgb(0,0,1,.2),border=NA)
  abline(v=input$cut_off[-1],lty=2,col="blue")
  abline(v=input$cut_off[1],col="blue")
  abline(v=lower,lty=2,col="green")
  abline(v=upper,col="green")
  polygon(c(lower, upper, upper, lower),c(0,0,1,1),col=rgb(0,1,0,.2),border=NA)
  dev.off()
}