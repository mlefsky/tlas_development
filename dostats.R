# go=function(){
#
#   code="656"
#
#   rm_nas=function(in1,in2){
#     #    print(in1)
#     print(in2)
#     fn
#     df <- data.frame(in1,in2)
#     names(df) = c("x","y")
#     tmp=(!is.na(df[,1])) & (!is.na(df[,2])) & (df[,1] != -9999) & ((df[,1] != -9999))
#     #    print(df)
#     #    print(tmp)
#     vars = subset(df, tmp)
#
#     #  ix1=(which((is.na(in1)==FALSE)|(in1=="NA")))
#     #  ix1=(which((is.na(in2)==FALSE)|(in2=="NA")))
#
#     return(vars)
#   }
#
#
#   if (code == "164"){tile_66_136_stack_edit=tile_66_136_164_stack}
#   if (code == "328"){tile_66_136_stack_edit=tile_66_136_328_stack}
#   if (code == "656"){tile_66_136_stack_edit=tile_66_136_656_stack}
#
#
#   tile_66_136_stack_edit[tile_66_136_stack_edit==-9999]= NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit==0]=NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit>1000000]= NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit<(-1000000)]= NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit==-Inf]=NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit==Inf]=NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit==Inf]=NA
#   tile_66_136_stack_edit[tile_66_136_stack_edit=="NA"]=NA
#
#   vnames=c("chm","cover","density_lt1","density",'pdsm',"pdtm")
#   dataset=c("sub_01","sub_02","sub_03","sub_04")
#   #print(dataset)
#   #zwl=list()
#
#   fname=paste("/home/lefsky/time_trials/trial_",code,"/table_data_",code,".csv",sep="")
#   df=data.frame()
#   #  with open('fname', 'w') as the_file:
#   
rm_nas = function(in1, in2) {
  df <- data.frame(in1, in2)
  names(df) = c("x", "y")
  tmp = (!is.na(df[, 1])) &
    (!is.na(df[, 2])) & (df[, 1] != -9999) & ((df[, 1] != -9999))
  #    print(df)
  #    print(tmp)
  vars = subset(df, tmp)
  
  #  ix1=(which((is.na(in1)==FALSE)|(in1=="NA")))
  #  ix1=(which((is.na(in2)==FALSE)|(in2=="NA")))
  
  return(vars)
}

go=function(code){
#  code="328"

if (code == "164") {
  tile_66_136_stack_edit = tile_66_136_164_stack
}
if (code == "328") {
  tile_66_136_stack_edit = tile_66_136_328_stack
}
if (code == "656") {
  tile_66_136_stack_edit = tile_66_136_656_stack
}


tile_66_136_stack_edit[tile_66_136_stack_edit == -9999] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit == 0] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit > 1000000] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit < (-1000000)] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit == -Inf] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit == Inf] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit == Inf] = NA
tile_66_136_stack_edit[tile_66_136_stack_edit == "NA"] = NA

vnames = c("chm", "cover", "density_lt1", "density", 'pdsm', "pdtm")
dataset = c("sub_01", "sub_02", "sub_03", "sub_04")
#print(dataset)
#zwl=list()

fname = paste("/home/lefsky/time_trials/trial_",
              code,
              "/table_data_",
              code,
              ".csv", 
              sep = "")
print(c("fn", fname))
df = data.frame()
#  with open('fname', 'w') as the_file:


for (vix in 1:6) {
  for (dix in 2:4) {
    vname1 = paste("sub_01", "_", vnames[[vix]], sep = "")
    vname2 = paste(dataset[dix], "_", vnames[[vix]], sep = "")#.replace("sub_sub","sub")
#    print(c(vname1, vname2))
#a    print(c(length(tile_66_136_stack_edit[[vname1]]),length(tile_66_136_stack_edit[[vname2]])))
    vars = rm_nas(tile_66_136_stack_edit[[vname1]], tile_66_136_stack_edit[[vname2]])
#    print(c(tile_66_136_stack_edit[[vname1]], tile_66_136_stack_edit[[vname2]]))
    #        #        print(vars)c()
    sd1 = sd(vars[["x"]] - vars[["y"]])
    x1 = mean(vars[["x"]] - vars[["y"]])
    xx = mean(vars[["x"]])
    xy = mean(vars[["y"]])
    line = paste(vname1, vname2, vix, dix, sd1, x1, xx, xy, sep = ",")
    print(line)
    #            the_file.write('Hello\n')     write.l(line, fileConn)
  }
}

fname = paste("/home/lefsky/time_trials/trial_",
              code,
              "/tile_66_136_",
              code,
              "_stack.csv",
              sep = "")
#print(c("fname", fname))
tmp = read.csv(fname)



if (code == "164") {
  data_164 = tmp
}
if (code == "328") {
  data_328 = tmp
}
if (code == "656") {
  data_656 = tmp
}}

go("164")
go("328")
go("656")


#   for (vix in 1:6) {
#         for (dix in 1:3) {
#           vname1=paste("sub_01","_",vnames[[vix]],sep="")
#           vname2=paste("sub_",dataset[dix],"_",vnames[[vix]],sep="")
#           print(c(vname1,vname2))
#           print(c(length(tile_66_136_stack_edit[[vname1]]),
#                   length(tile_66_136_stack_edit[[vname2]])))
#           vars=rm_nas(tile_66_136_stack_edit[[vname1]],tile_66_136_stack_edit[[vname2]])
#           #        print(vars)
#           sd1=sd(vars[["x"]]-vars[["y"]])
#           x1=mean(vars[["x"]]-vars[["y"]])
#           xx=mean(vars[["x"]])
#           xy=mean(vars[["y"]])
#           line=paste(vname1,vname2,vix,dix,sd1,x1,xx,xy,sep=",")
#           print(line)
# #            the_file.write('Hello\n')     write.l(line, fileConn)
#           }
#     }
#
#     close(fileConn)
#
# fname = paste("/home/lefsky/time_trials/trial_",
#               code,
#               "/tile_66_136_",
#               code,
#               "_stack.csv",
#               sep = "")
# print(c("fname2: ", fname))
# tmp = read.csv(fname)
# 
# 
# 
# if (code == "164") {
#   data_164 = tmp
# }
# if (code == "328") {
#   data_328 = tmp
# }
# if (code == "656") {
#   data_656 = tmp
# }
# 


