my_ggplot1 = function(df,a,mean0,mean1,mean2,sd0,sd1,sd2,Amin,
                      lsize=12,lunit=0.8,tsize=11)
{
  lab = c(paste0(sprintf("%.2f",mean0[a]),"(",
               sprintf("%.2f",sd0[a]),")"),
          paste0(sprintf("%.2f",mean1[a]),"(",
               sprintf("%.2f",sd1[a]),")"),
          paste0(sprintf("%.2f",mean2[a]),"(",
               sprintf("%.2f",sd2[a]),")"))
  ggplot(df, aes(x = y, color = Method, fill = Method)) +
  geom_density(alpha = 0.5) +
  ggtitle(paste0("Age = ", a + Amin - 1)) +
  xlab("") + ylab("") +
  # Remove legend title + set custom labels
  scale_color_discrete(
    name = NULL,
    labels = lab) +
  scale_fill_discrete(
    name = NULL,
    labels = lab)+
  # Remove ANY remaining legend title (merged guide)
  guides(
    color = guide_legend(title = NULL),
    fill  = guide_legend(title = NULL)) +
  theme(legend.position = c(0.98,0.98),legend.justification = c(1, 1),
        legend.text=element_text(size=lsize),
        legend.key.size = unit(lunit, "lines"),
        legend.title = element_blank(),
        plot.margin = margin(20, 20, 20, 20),
        legend.box.margin = margin(10, 10, 10, 10),
        plot.title=element_text(size=tsize))
}

my_ggplot2 = function(df,a,mean0,mean1,msd0,sd1,Amin,
                      lsize=12,lunit=0.8,tsize=11)
{
  lab = c(paste0(sprintf("%.2f",mean0[a]),"(",
                 sprintf("%.2f",sd0[a]),")"),
          paste0(sprintf("%.2f",mean1[a]),"(",
                 sprintf("%.2f",sd1[a]),")"))
  ggplot(df, aes(x = y, color = Method, fill = Method)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste0("Age = ", a + Amin - 1)) +
    xlab("") + ylab("") +
    # Remove legend title + set custom labels
    scale_color_discrete(
      name = NULL,
      labels = lab) +
    scale_fill_discrete(
      name = NULL,
      labels = lab )+
    # Remove ANY remaining legend title (merged guide)
    guides(
      color = guide_legend(title = NULL),
      fill  = guide_legend(title = NULL)) +
    theme(legend.position = c(0.98,0.98),legend.justification = c(1, 1),
          legend.text=element_text(size=lsize),
          legend.key.size = unit(lunit, "lines"),
          legend.title = element_blank(),
          plot.margin = margin(20, 20, 20, 20),
          legend.box.margin = margin(10, 10, 10, 10),
          plot.title=element_text(size=tsize))
}