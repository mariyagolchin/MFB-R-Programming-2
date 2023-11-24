

PlotHeatmap<-function (matrix_gene_sample,gene_nemes,matrix_gene_variants){
  Ylength <- ncol(matrix_gene_sample) #193
  print(Ylength)
  Xlength <- nrow(matrix_gene_sample) #3 different with number of rows that enter
  # maxColorSize = 14
  plot (c(-10,Ylength+10),c(-10,Xlength+20),
        type= 'n',
        axes = FALSE,
        ann=FALSE)
  variations_name = c(
  'Frame_Shift_Del'=1,
  'Missense_Mutation'=2,
  'Nonsense_Mutation'=3,
  'Frame_Shift_Ins'=4,
  'In_Frame_Ins'=5,
  'Splice_Site'=6,
  'In_Frame_Del'=7,
  'Multi_Hit'=8
  )
  colors = c('gray','blue','green','brown','purple','red','orange','yellow','black')
  for(i in Xlength:1){
        for (j in 1:Ylength){
          x= matrix_gene_sample[i,j]
          print(paste("x=>",x,"col=>",colors[x+1]))
          xCoords = c(j,j+1,j+1,j)
          yCoords = c(i,i,i-1,i-1)
          polygon(xCoords,yCoords,col = colors[x+1])
        }
        text(-10,i-0.75,gene_nemes[i],cex = 0.5)
        # print("=================next gene================matrix_gene_variants rows 193 samples ,col 8 variation===============")
        # drawRightBoxPlot(i,j+10,rightBar[i-1,],max(matrix_gene_variants))
        variations= matrix_gene_variants[i,]
        print(matrix_gene_variants[i,])
        for (variat in variations){
          print(variat)
          colors_part_each_variation =(variat/ sum(variations)*10)
          print(paste("colors_part_each_variation==>",colors_part_each_variation))
          xCoords = c(j+5, j+5+colors_part_each_variation, j+5+colors_part_each_variation, j+5)
          yCoords = c(i-0.1,i-0.1,i-1+0.1,i-1+0.1)
          co = c('blue','green','brown','purple','red','orange','yellow','black')
          polygon(xCoords,yCoords,col = co[variat])
          j=j+colors_part_each_variation
        }
    
  }
}
###################read data######################################
##################################################################
data_new= read.csv("tcga_laml.maf",sep = '\t') #dim=> 2207   17
##### call function matrix_gene_sample
View(data_new)
# ========== matrix_gene_sample ===================
# Gene_Sample <- function(data_new){
matrix_gene_sample<- matrix(0,length(unique(data_new$Hugo_Symbol)), length(unique(data_new$Tumor_Sample_Barcode)))
rownames(matrix_gene_sample) = unique(data_new$Hugo_Symbol)
colnames(matrix_gene_sample) = unique(data_new$Tumor_Sample_Barcode)
variants<-unique(data_new$Variant_Classification)
variations_name = c(
  'Frame_Shift_Del'=1,
  'Missense_Mutation'=2,
  'Nonsense_Mutation'=3,
  'Frame_Shift_Ins'=4,
  'In_Frame_Ins'=5,
  'Splice_Site'=6,
  'In_Frame_Del'=7,
  'Multi_Hit'=8
)
for(i in 1:dim(data_new)[1] ){
  print(i)  
  hugo_i<- data_new$Hugo_Symbol[i]
  sample_i <- data_new$Tumor_Sample_Barcode[i]
  variant_i<- data_new$Variant_Classification[i]
  print(paste(hugo_i,",", sample_i,", ",variant_i))
  if (variant_i %in% names(variations_name)){
    print(variations_name[[variant_i]])
    if (matrix_gene_sample[hugo_i, sample_i] == 0){
      matrix_gene_sample[hugo_i, sample_i] = variations_name[[variant_i]]
    }else{
      if (matrix_gene_sample[hugo_i, sample_i] != variations_name[[variant_i]])
      {
        matrix_gene_sample[hugo_i, sample_i] = variations_name[['Multi_Hit']]
      }
    }
  }
}

######### sorted_matrix_gene_sample by_mostvariants ##### sort matrix by gene with most variant at tops  ##############
sorted_matrix_by_mostvariants<- matrix_gene_sample[order(apply(matrix_gene_sample, 1, function(row) sum(row != 0)), decreasing = T), ] # 
sum(sorted_matrix_by_mostvariants[1,] !=0)# 52
sum(sorted_matrix_by_mostvariants[2,] !=0) # 48
############right plot########### matrix_gene_variants########################################

matrix_gene_variants = matrix(0,length(unique(data_new$Hugo_Symbol)),length(variations_name) )
dim(matrix_gene_variants) #1611   8
# dim(sorted_matrix_by_mostvariants)#1611  193
for (i in 1:dim(sorted_matrix_by_mostvariants)[1])#1611 gene
{
  for (j in 1: dim(sorted_matrix_by_mostvariants)[2])#193 sample
  {
    variant_i_j = sorted_matrix_by_mostvariants[i,j] # variant gene_i, sample_j
    matrix_gene_variants[i,variant_i_j] = matrix_gene_variants[i,variant_i_j] + 1
  }
}

# ===check===
sorted_matrix_by_mostvariants[1,]
sum(sorted_matrix_by_mostvariants[1,] !=0)#52
matrix_gene_variants[1,]
# 0 15  0  0 33  3  1  0
# draw barplot
plot(c(-3,50),c(-3,50))
for ( item in matrix_gene_variants[1,]){
    counter= 1
    print(item)
    polygon(c(counter, counter, counter+0.05, counter+0.05), c(0, item, item, 0))
    text(x-1, -3, item , cex=0.5, srt=90)
    counter=counter+2

}

sum(matrix_gene_variants[1,]) #52
##########################################################################################
################ sort and plot with no colere ###################

GS<-sorted_matrix_by_mostvariants[1:30,]
gs_01<-GS[which(GS == 1)] <- 0
gs_01<-GS[which(GS != 0)] <- 1
# T<-order(GS[1,],GS[2,],GS[3,],GS[4,],GS[5,],decreasing=TRUE)

T<-order(GS[1,],GS[2,],GS[3,],GS[4,],GS[5,],GS[6,],GS[7,],GS[8,],GS[9,],GS[10,],
         GS[11,],GS[12,],GS[13,],GS[14,],GS[15,],GS[16,],GS[17,],GS[18,],GS[19,],GS[20,]
         ,GS[21,],GS[22,],GS[23,],GS[24,],GS[25,],GS[26,],GS[27,],GS[28,],GS[29,],GS[30,],decreasing=TRUE)

GS[,T]
gene_nemes=rownames(sorted_matrix_by_mostvariants)[1:30]

# pdf('PlotHeatmap_nocolor.pdf')
PlotHeatmap(GS[,T],gene_nemes,matrix_gene_variants)
# dev.off()

# ================= sort and plot with colere===============================
GS<-sorted_matrix_by_mostvariants[1:30,]
# gs_01<-GS[which(GS == 1)] <- 0
# gs_01<-GS[which(GS != 0)] <- 1
# GS
# T<-order(GS[1,],GS[2,],GS[3,],GS[4,],GS[5,],decreasing=TRUE)

T<-order(GS[1,],GS[2,],GS[3,],GS[4,],GS[5,],GS[6,],GS[7,],GS[8,],GS[9,],GS[10,],
         GS[11,],GS[12,],GS[13,],GS[14,],GS[15,],GS[16,],GS[17,],GS[18,],GS[19,],GS[20,]
         ,GS[21,],GS[22,],GS[23,],GS[24,],GS[25,],GS[26,],GS[27,],GS[28,],GS[29,],GS[30,],decreasing=TRUE)

GS[,T]

# pdf('PlotHeatmap_with_color.pdf')
gene_nemes=rownames(sorted_matrix_by_mostvariants)[1:30]
PlotHeatmap(GS[,T],gene_nemes, matrix_gene_variants)
# dev.off()


# ====== top plot ==== matrix_sample_variants ===================
# matrix_sample_variants = matrix(0,length(unique(data_new$Tumor_Sample_Barcode)),length(variations_name)-1)
# dim(matrix_sample_variants) # 193  7
# for(i in 1:dim(data_new)[1] ){
#   print(i)
#   hugo_i<- data_new$Hugo_Symbol[i]
#   sample_i <- data_new$Tumor_Sample_Barcode[i]
#   variant_i<- data_new$Variant_Classification[i]
#   print(paste(hugo_i,",", sample_i,", ",variant_i))
#   if (variant_i %in% names(variations_name)){
#     matrix_sample_variants[i,variant_i_j] = matrix_sample_variants[i,variant_i_j] + 1
#   }
#   }
