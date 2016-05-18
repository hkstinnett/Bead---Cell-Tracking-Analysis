#save project folder
old_wd <- getwd()

# This isn't a command, just reminder to change the working directory
# to a data folder so the data can be loaded up.
Setwd()



ID<-dir()

# fix header issue from Fiji in data files, renames columns

for (i in 1:(length(ID)/2)){
  tmpbead<-readLines(ID[2*i-1])
  tmpbead[1]<-"Serial\tTrack\tSlice\tX\tY\tDistance\tVelocity\tPixel"
  writeLines(tmpbead,con=ID[2*i-1])
  tmpcells<-readLines(ID[2*i])
  tmpcells[1]<-"Serial\tTrack\tSlice\tX\tY\tDistance\tVelocity\tPixel"
  writeLines(tmpcells,con=ID[2*i])
}



#make lists to sort data into

bead.fgf8<-list(NULL)
cells.fgf8<-list(NULL)
bead.ctrl<-list(NULL)
cells.ctrl<-list(NULL)
j=0
k=0

#import data files and group in the lists according to name

for (i in 1:(length(ID)/2)){
  if(charmatch("Ctrl", substr(ID[2*i-1],1,4),nomatch=0)){
    j=j+1
    bead.ctrl[[j]]<-read.table(file=ID[2*i-1],head=T)
    cells.ctrl[[j]]<-read.table(file=ID[2*i],head=T)
    print(paste("ctrl:",ID[2*i-1],ID[2*i]))
  }
  else{
    k=k+1
    bead.fgf8[[k]]<-read.table(file=ID[2*i-1],head=T)
    cells.fgf8[[k]]<-read.table(file=ID[2*i],head=T)
    print(paste("fgf8:",ID[2*i-1],ID[2*i]))
  }
}

# head back to project folder
setwd(old_wd)
#make lists for subtracted tracks

sub.fgf8<-list(NULL)
sub.ctrl<-list(NULL)

#subtracting bead distance from cell distance for each time point per cell
#unnecessary to use rep, R will continue until the sets are even
#yeah

for (i in 1:(length(cells.fgf8))){
  sub.fgf8[[i]]<-data.frame(X=(cells.fgf8[[i]]$X-bead.fgf8[[i]]$X),
                            Y=(cells.fgf8[[i]]$Y-bead.fgf8[[i]]$Y),
                            Time=cells.fgf8[[i]]$Slice,
                            Track=cells.fgf8[[i]]$Track)
}

for (i in 1:(length(cells.ctrl))){
  sub.ctrl[[i]]<-data.frame(X=(cells.ctrl[[i]]$X-bead.ctrl[[i]]$X),
                            Y=(cells.ctrl[[i]]$Y-bead.ctrl[[i]]$Y),
                            Time=cells.ctrl[[i]]$Slice,
                            Track=cells.ctrl[[i]]$Track)
}

#this number 30 is hard coded, should be changed to 20 in fgf or we should figure out 
#how to fix this one
#commenting out below to show what was changed from the previous version. Don't need to use rep.

#for (i in 1:(length(cells.fgf8))){
 # sub.fgf8[[i]]<-data.frame(X=(cells.fgf8[[i]]$X-rep(bead.fgf8[[i]]$X,20)),
  #                          Y=(cells.fgf8[[i]]$Y-rep(bead.fgf8[[i]]$Y,20)),
   #                         Time=cells.fgf8[[i]]$Slice,
    #                        Track=cells.fgf8[[i]]$Track)
#}

#for (i in 1:(length(cells.ctrl))){
 # sub.ctrl[[i]]<-data.frame(X=(cells.ctrl[[i]]$X-rep(bead.ctrl[[i]]$X,20)),
  #                          Y=(cells.ctrl[[i]]$Y-rep(bead.ctrl[[i]]$Y,20)),
   #                         Time=cells.ctrl[[i]]$Slice,
    #                        Track=cells.ctrl[[i]]$Track)
#}

#set n timepoints equal to length of tracks in bead, always 100
#divide cell serials by n_timepoints to calculate # of cells
#assumes all experiments have same number of timepoints

n_fgftimepoints <- length(bead.fgf8[[1]]$Track)
n_fgfcells <- length(cells.fgf8[[1]]$Serial)/n_fgftimepoints


Distance.fgf8<-sapply(sub.fgf8,
                      function (embryo) by(embryo, embryo$Track, 
                      function(embryo.track) sqrt(embryo.track$X^2+embryo.track$Y^2)))
Distance.fgf8<-unlist(Distance.fgf8)
dim(Distance.fgf8)<-c(n_fgftimepoints, n_fgfcells,(length(sub.fgf8)))

n_ctrltimepoints <- length(bead.ctrl[[1]]$Track)
n_ctrlcells <- length(cells.ctrl[[1]]$Serial)/n_ctrltimepoints

Distance.ctrl<-sapply(sub.ctrl, 
                      function (embryo) by(embryo, embryo$Track, 
                      function(embryo.track) sqrt(embryo.track$X^2+embryo.track$Y^2)))
Distance.ctrl<-unlist(Distance.ctrl)
dim(Distance.ctrl)<-c(n_ctrltimepoints,n_ctrlcells,(length(sub.ctrl)))

#fgf with the # 20

#Distance.fgf8<-sapply(sub.fgf8,
 #                     function (embryo) by(embryo, embryo$Track, 
 #                     function(embryo.track) sqrt(embryo.track$X^2+embryo.track$Y^2)))
#Distance.fgf8<-unlist(Distance.fgf8)
#dim(Distance.fgf8)<-c(100,20,(length(sub.fgf8)))

#Distance.ctrl<-sapply(sub.ctrl, 
 #                     function (embryo) by(embryo, embryo$Track, 
  #                    function(embryo.track) sqrt(embryo.track$X^2+embryo.track$Y^2)))
#Distance.ctrl<-unlist(Distance.ctrl)
#dim(Distance.ctrl)<-c(100,20,(length(sub.ctrl)))



Distance.fgf8.mean<-apply(Distance.fgf8, 1, mean)
Distance.ctrl.mean<-apply(Distance.ctrl, 1, mean)
Distance.fgf8.sd<-apply(Distance.fgf8, 1, sd)
Distance.ctrl.sd<-apply(Distance.ctrl, 1, sd)
par(mfrow=c(1,2))

#to plot distances in microns

plot(1:100, Distance.fgf8.mean-Distance.fgf8.mean[1],
     type="l", lwd=2, col="maroon",main="Fgf24 deficient background",
     xlab="Time/8 min",ylab="change in average bead-cell distance (um)", ylim=c(-50, 10))


lines(1:100, lwd=2, Distance.ctrl.mean-Distance.ctrl.mean[1], col="black")



#180 is the number of points being averaged over "n" calculated from shape of Distance.fgf8.mean dim
fgf8_N<-prod(dim(Distance.fgf8)[2:3])

CI_high_fgf<-(Distance.fgf8.mean+1.96*Distance.fgf8.sd/sqrt(fgf8_N))-Distance.fgf8.mean[1]
CI_low_fgf<-(Distance.fgf8.mean-1.96*Distance.fgf8.sd/sqrt(fgf8_N))-Distance.fgf8.mean[1]
polygon(c(1:100, rev(1:100)), c(CI_high_fgf, rev(CI_low_fgf)), col = rgb(.6, 0, .2, alpha=0.2), border = NA)


ctrl_N<-prod(dim(Distance.ctrl)[2:3])

CI_high_ctrl<-(Distance.ctrl.mean+1.96*Distance.ctrl.sd/sqrt(ctrl_N))-Distance.ctrl.mean[1]
CI_low_ctrl<-(Distance.ctrl.mean-1.96*Distance.ctrl.sd/sqrt(ctrl_N))-Distance.ctrl.mean[1]
polygon(c(1:100, rev(1:100)), c(CI_high_ctrl, rev(CI_low_ctrl)), col = rgb(.2, .2, .2, alpha=0.2), border = NA)


#to plot fold changes

#commenting this out, so I can use the ylim specs below. ylim to change axes
#plot(1:100, Distance.fgf8.mean/Distance.fgf8.mean[1],
#    type="l",col="maroon",main="fgf8 in maroon, control in black",
#   xlab="Time/8 min",ylab="fold change of bead-cell distance")

plot(1:100, Distance.fgf8.mean/Distance.fgf8.mean[1],type="l",col="maroon",main="fgf8 in maroon, control in black",xlab="Time/8 min",ylab="fold change of bead-cell distance",
ylim=c(0.5,1.5))
#use ylim= if you need to change axes
lines(1:100, Distance.ctrl.mean/Distance.ctrl.mean[1], col="black")
wilcox.test(Distance.ctrl,Distance.fgf8)

CI_high_fgf<-(Distance.fgf8.mean+1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1]
CI_low_fgf<-(Distance.fgf8.mean-1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1]

polygon(c(1:100, rev(1:100)), c(CI_high_fgf, rev(CI_low_fgf)), col = rgb(0.7, 0, 0.2, alpha=0.1), border = NA)

CI_high_ctrl<-(Distance.ctrl.mean+1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1]
CI_low_ctrl<-(Distance.ctrl.mean-1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1]

polygon(c(1:100, rev(1:100)), c(CI_high_ctrl, rev(CI_low_ctrl)), col = rgb(0.2, 0.2, 0.2, alpha=0.1), border = NA)

#error lines with dashes
#lines(1:100, (Distance.fgf8.mean+1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1], col="maroon",lty="dashed")
#lines(1:100, (Distance.fgf8.mean-1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1], col="maroon",lty="dashed")
#lines(1:100, (Distance.ctrl.mean+1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1], col="black",lty="dashed")
#lines(1:100, (Distance.ctrl.mean-1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1], col="black",lty="dashed")

#to plot a non-normalized version and see the 
#distances of the mean from the bead over time

plot(1:100, Distance.fgf8.mean, type="l",col="maroon",
     main="fgf8 in maroon, control in black", 
     xlab="Time/8 min",ylab="average bead-cell distance (?m)", 
     ylim = c(72, 130) )
lines(1:100, Distance.ctrl.mean, col="black")

#I want to plot each individual distance, not just mean.

plot(1:100, Distance.fgf8[1:5], type="l", col="maroon",xlab="Time/8 min",ylab="average bead-cell distance (?m)")


#lines(1:100, (Distance.fgf8.mean+1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1], col="orange",lty="dashed")
#lines(1:100, (Distance.fgf8.mean-1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1], col="orange",lty="dashed")

#add CI for this data


CI_high_fgf<-(Distance.fgf8.mean+1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1]
CI_low_fgf<-(Distance.fgf8.mean-1.96*Distance.fgf8.sd/sqrt(160))/Distance.fgf8.mean[1]
polygon(c(1:100, rev(1:100)), c(CI_high_fgf, rev(CI_low_fgf)), col = rgb(.6, 0, .2, alpha=0.2), border = NA)
lines(1:100, Distance.fgf8.mean/Distance.fgf8.mean[1], col="maroon")

CI_high_ctrl<-(Distance.ctrl.mean+1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1]
CI_low_ctrl<-(Distance.ctrl.mean-1.96*Distance.ctrl.sd/sqrt(160))/Distance.ctrl.mean[1]
polygon(c(1:100, rev(1:100)), c(CI_high_ctrl, rev(CI_low_ctrl)), col = rgb(0.2,0.2,0.2, alpha=0.2), border = NA)
lines(1:100, Distance.ctrl.mean/Distance.ctrl.mean[1], col="black")


