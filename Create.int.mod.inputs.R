fn.paste.vec=function(txt) paste(txt, collapse='\t')
fn.paste.c=function(val,txt) paste(val, paste("\t\t#",txt))
fn.col.nms.c=function(NMS) c('\n#',paste(NMS,collapse='\t'),'\n')
create.ass.inputs=function(WD,Outfile,Input_list,extension=".DAT")
{
  tmp=list()
  for(u in 1:length(Input_list))
  {
    dummy=any(is.data.frame(Input_list[[u]]),is.matrix(Input_list[[u]]))
    
    if(dummy)
    {
      for(d in 1:nrow(Input_list[[u]])) tmp <- c(tmp, paste(Input_list[[u]][d,], collapse = "\t"),"\n")
    }else
    {
      tmp=c(tmp,Input_list[[u]])
    }
  }
    
  write.table(tmp, paste(WD,paste(Outfile,extension,sep=''),sep='/'),
              sep="", row.names = F, col.names = F, quote=F)
}