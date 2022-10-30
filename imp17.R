library("RColorBrewer")
library(tidyverse)
library(ggplot2)
library(GGally)
library(grid)
library(gridExtra)

p_<-GGally::print_if_interactive
points_legend<-gglegend(ggally_points)

#loading the dataset
load("/cloud/project/df-forAranyak.RData")


print(df[,6])

n<-nrow(df)

X_hat<-cbind(df$X.1,df$X.2,df$X.3,df$X.4,df$X.5,df$X.6)
y<-df$dist
claw<-as.character(df$claw)

ggpairs(df[,6:11],
        lower = list(continuous="points"),
        upper = list(continuous=gglegend("points")),
        mapping=ggplot2::
          aes(colour = claw))+
  scale_x_continuous(guide = 
                       guide_axis(n.dodge = 1)) +
  scale_y_continuous(guide = 
                       guide_axis(check.overlap = 
                                    TRUE)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_colour_manual(values = c("red","purple",
                                 "yellow","pink",
                                 "green","brown",
                                 "blue"),
                      labels=unname(TeX(c(
                        "claw-0",
                        "claw-1",
                        "claw-2",
                        "claw-3",
                        "claw-4",
                        "claw-5",
                        "claw-6")))) +
  legend("topright",legend = c(
    "claw-0",
    "claw-1",
    "claw-2",
    "claw-3",
    "claw-4",
    "claw-5",
    "claw-6"),
    col = c("red","purple",
            "yellow","pink",
            "green","brown",
            "blue"))

points_legend <- gglegend(ggally_points)
p_(points_legend(df[,6:11],ggplot2::aes(df[,6], 
                                    df[,7], 
                                    color = 
                                      as.character(df$claw)
                                    )))

pm<-ggpairs(df[,6:11],
            mapping = ggplot2::aes(color = as.character(df$claw)),
            upper = list(continuous = gglegend("points"),
            )


cols<-character(nrow(df))
cols[]<-"black"
print(cols)

cols[df$claw==0]<-"red"
cols[df$claw==1]<-"orange"
cols[df$claw==2]<-"blue"
cols[df$claw==3]<-"green"
cols[df$claw==4]<-"yellow"
cols[df$claw==5]<-"gray"
cols[df$claw==6]<-"purple"
cols[df$claw==7]<-"brown"

pairs(df[,6:11],col=cols)