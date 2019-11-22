library(ggplot2)
library(forcats)
library(stringr)


A<-read.table(file="/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/GSEA_cell_data_new/deci_T_hmk_new/gsea_report_for_deciNM_1573569795364.xls",sep="\t",head=T,row.names=1,fill=T)
B<-read.table(file="/home/yjingjing/project/hongfangzi/hongfangzi_MATRIXrawdata/GSEA_cell_data_new/deci_T_hmk_new/gsea_report_for_deciPE_1573569795364.xls",sep="\t",head=T,row.names=1,fill=T)
gsea_dot <- rbind(A,B)
dim(gsea_dot)
gsea_dot <- gsea_dot[which(gsea_dot$FDR.q.val < 0.05), ]
dim(gsea_dot)
p<- ggplot(gsea_dot, aes(x = NES, y = fct_reorder(GS.br..follow.link.to.MSigDB, NES))) + 
  geom_point(aes(color = FDR.q.val,size = SIZE)) +
  coord_cartesian(xlim=c(-3,3))+ 
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3))+
  #
  scale_colour_gradientn(limits=c(0, 1), colours=rainbow(6)) +
  theme_bw(base_size = 10) +
  ylab(NULL) +scale_fill_brewer(palette = 'Accent')+theme(panel.background=element_rect(fill='grey95'),panel.border = element_blank())+
  geom_vline(xintercept = 0,size=1)+
  geom_segment(aes(x=-2,y=0,xend=2,yend=0)) +

  theme(axis.text.y=element_text(size=15),axis.text.x = element_text(size = 7))+
  scale_y_discrete(labels=function(x)str_wrap(x,width=35))   

pdf(file="./dotplot1_gsea_deci_T.pdf",8,8) 
p
dev.off() 
