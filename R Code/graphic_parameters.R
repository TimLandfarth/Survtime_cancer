# 1. Graphic parameters----
## 1.1 Colors----
gray1 <- "#969696"
gray2 <- c("#636363","#969696")
gray3 <- c("#636363","#969696","#bdbdbd")
col_s <- c("#8C8C8C", "#88BDE6", "#FBB258", "#90CD97", "#F6AAC9", "#BFA554", "#A899C7", "#EDDD46", "#F07E6E")
col_m <- c("#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF35", "#F15854")
col_h <- c("#000000", "#265DAB", "#DF5C24", "#059748", "#E5126F", "#9D722A", "#7B3A96", "#C7B42D", "#CB2027")
axis_col <- "grey70"

## 1.1 Sizes----
linesize <- 1.5
pointsize <- 5
axissize <- .5
pdf_w_h <- c((7/3)*4, (7/3) * 3)

# 2. functions----
## 2.1 Outcome vs. covariate----
plot_outcov <- function(data = df, cov, covname) {
  require(ggplot2)
  require(ggthemes)
  require(rlang)
  return(
    ggplot(data, aes(y =  .data[[cov]] , x = time)) +
      geom_point() +
      theme_tufte() +
      theme(
        axis.line.x = element_line(gray2[1]),
        axis.line.y = element_line(gray2[1])
      ) +
      labs(y = covname,
           x = "time")
  )
}
