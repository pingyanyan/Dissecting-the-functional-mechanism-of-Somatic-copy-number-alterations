

common_data<-function(cnv_M,gene_M,lnc_M,mir_M){

    common_sample<-intersect(intersect(colnames(cnv_M),colnames(gene_M)),intersect(colnames(lnc_M),colnames(mir_M)))
    
	cnv_M<-cnv_M[,common_sample]
	gene_M<-gene_M[,common_sample]
	lnc_M<-lnc_M[,common_sample]
	mir_M<-mir_M[,common_sample]
	
	common_gene<-c(intersect(rownames(cnv_M),rownames(gene_M)),intersect(rownames(cnv_M),rownames(lnc_M)))
	cnv_M<-cnv_M[common_gene,]##
	save(cnv_M,gene_M,lnc_M,mir_M,file="common_data.Rdata")
}

######################################################################################################################
##
del_cnv<-function(cnv_M, freq_thr){
      pos<-which(rowSums(cnv_M!=0)>=dim(cnv_M)[2]*freq_thr)
	  cnv_M<-cnv_M[pos,]
	  return(cnv_M)
	  
}

######################################################################################################################
##3.

del_exp<-function(M,exp_thr,freq_thr){
  pos1<-which(apply(M,1,mean)>=exp_thr)
  pos2<-which(rowSums(M==0)<=dim(M)[2]*freq_thr)
  pos<-intersect(pos1,pos2)
  M<-M[pos,]
  return(M)
}
##########################################################################################################################

#
######################################################################################################################
#####

common_data1<-function(M,inter_list){

	 common_gene<-intersect(rownames(M),names(inter_list))
	 M<-M[common_gene,]
	 if(length(common_gene)){
	    return(M)
	 }else{
	    return(NULL)
	 }

}
######################################################################################################################
#####
common_data2<-function(M,inter_list){
 
	 inter_list<-lapply(inter_list,intersect,y=rownames(M))
	 l<-unlist(lapply(inter_list,length))
	 pos<-which(l==0)
	 if(length(pos)){
	    inter_list<-inter_list[-pos]
	 }
	 return(inter_list)

}

######################################################################################################################################################
###
cnvs_exp<-function(cnv_M,gene_exp,p_thr,min_thr){
 
     common_gene<-as.matrix(intersect(rownames(cnv_M),rownames(gene_exp)))
	 result<-apply(common_gene,1,cnv_exp2,cnv_M=cnv_M,gene_exp=gene_exp,p_thr=p_thr,min_thr=min_thr)
	
	 pos<-which(result)
	 result<-common_gene[pos,1]
	 return(result)
}
 
cnv_exp1<-function(cnv_flag,gene_exp,p_thr,min_thr){

    sample_tag<-tapply(cnv_flag,as.factor(cnv_flag),names)
	pos<-which(unlist(lapply(sample_tag,length))<min_thr) 
	if(length(pos)){
	  sample_tag<-sample_tag[-pos]
	}
	
	if(length(sample_tag)==2){ 
	     tt<-t.test(gene_exp[sample_tag[[1]]],gene_exp[sample_tag[[2]]]) 
	     p<-tt$p.value
		 FD<-(tt$estimate)[2]-(tt$estimate)[1]
		 tag<-(p<=p_thr&FD>0)
	
	}else{ 
	     tt1<- t.test(gene_exp[sample_tag[[1]]],gene_exp[sample_tag[[2]]])
		 
		 p<-tt1$p.value
		 FD<-(tt1$estimate)[2]-(tt1$estimate)[1]
		 
		 tag1<-(p<=p_thr&FD>0)
		 
	     tt2<-t.test(gene_exp[sample_tag[[2]]],gene_exp[sample_tag[[3]]])
	     p<-tt2$p.value
		 FD<-(tt2$estimate)[2]-(tt2$estimate)[1]
		 
		 tag2<-(p<=p_thr&FD>0)
		 tag<-(tag1|tag1)
	
	}
	return(tag) 

}
##############################
###
cnv_exp2<-function(geneid,cnv_M,gene_exp,p_thr,min_thr){

      result<-cnv_exp1(cnv_M[geneid,],gene_exp[geneid,],p_thr,min_thr) 
	  return(result)
}
#############################################################################################################################
##############################
###
ceRNA_mir_ceRNA<-function(ceRNA,ceRNA_list,mir_list){
      
	   mirs<-ceRNA_list[[ceRNA]]
	   common_mirs<-intersect(mirs,names(mir_list))
	   temp_list<-mir_list[common_mirs]
	   lens<-unlist(lapply(temp_list,length))
	   mir_ceRNA<-cbind(rep(names(temp_list),lens),unlist(temp_list))
	   l<-dim(mir_ceRNA)[1]
	   ceRNA_mir_ceRNAs<-cbind(rep(ceRNA,l),mir_ceRNA)
	   
	   return(ceRNA_mir_ceRNAs)
}
###########################################################################
######################################################################################################################
##
profile_class<-function(flag,gene_exp,mir_M){
## 
       a_list<-list()
		length(a_list)<-2
		names(a_list)<-c("gene","mir")
		
		a_list$gene<-gene_exp[,flag]
		a_list$mir<-mir_M[,flag]
		return(a_list)
		

}
################################################################
##
case_normal_profile<-function(flag_v,gene_exp,mir_M,min_thr){
 
 
		alter_sample<-tapply(flag_v,as.factor(flag_v),names)
		
		pos<-which(unlist(lapply(alter_sample,length))<min_thr)
		if(length(pos)){
		  alter_sample<-alter_sample[-pos]
		 }
		 
        result<-lapply(alter_sample,profile_class,gene_exp=gene_exp,mir_M=mir_M)
		return(result)
          
}


#######################################################################################
##
cors<-function(a_pair,expr1,expr2){
 
       temp<-cor.test(expr1[a_pair[1],],expr2[a_pair[2],])
	   pcor<-temp$estimate
	   p<-temp$ p.value
	   result<-c(pcor,p)
	   return(result)

}

corrs<-function(all_pairs,expr1,expr2){
## 
    c_p<-t(apply(all_pairs,1,cors,expr1=expr1,expr2=expr2))
    colnames(c_p)<-c("cor","pvalue")
	c_p<-cbind(all_pairs,as.data.frame(c_p))
	return(c_p)
}
########################################################################################

### 
pair_filter<-function(c_p,flag,method,cor_thr,pthr){
 ### 
    if(flag=="mir"){ 
	   pos1<-which(c_p[,3]<=cor_thr)
	}else{
	   pos1<-which(c_p[,3]>=cor_thr)
	}

    Padj<-p.adjust(c_p[,4],method=method)
	c_p<-cbind(c_p,Padj)
	Padj<-c_p[,5]
	pos2<-which(Padj<=pthr) 
	pos<-intersect(pos1,pos2)
	return(c_p[pos,])
}
########################################################################

 
filter1<-function(ceRNA,sig_pair,flag){
### 
     if(length(ceRNA)==0|length(sig_pair)==0){
	    return(NULL)
	 }else if(length(ceRNA)==3){
	     ceRNA<-matrix(ceRNA,nrow=1)
	 
	 }
        sig_pairs<-rownames(sig_pair)
        if(flag=="12"){
		  temp_pair<-paste(ceRNA[,1],ceRNA[,2],sep="_")
		  
		}else if(flag=="13"){
		  temp_pair<-paste(ceRNA[,1],ceRNA[,3],sep="_")
		}else{
		  temp_pair<-paste(ceRNA[,2],ceRNA[,3],sep="_")
		}
		pos<-which(temp_pair%in%sig_pairs)
		
		return(ceRNA[pos,])

}

#############################################################

### 
extract_cp<-function(ceRNA,flag,c_ps){
 
    if(length(ceRNA)==0){
	    return(NULL)
	 }else if(length(ceRNA)==3){
	     ceRNA<-matrix(ceRNA,nrow=1)
	 
	 }



        sig_pairs<-rownames(c_ps)
		
        if(flag=="12"){
           name<-paste(ceRNA[,1],ceRNA[,2],sep="_")
		}else if(flag=="13"){
		    name<-paste(ceRNA[,1],ceRNA[,3],sep="_")
		}else{
		   name<-paste(ceRNA[,2],ceRNA[,3],sep="_")
		}
        #result<-cbind(ceRNA,c_ps[name,c("cor","pvalue","Padj")])
		result<-cbind(ceRNA,c_ps[name,c("cor","pvalue")])
		return(result)
}
#################################################################

cor_all<-function(ceRNA_triples,mir_M,gene_exp,cthr1,cthr2,cthr3,pthr){
####### 
  mir_ceRNA1<-unique(ceRNA_triples[,1:2])
  colnames(mir_ceRNA1)<-c("ceRNA1","mir")
  mir_ceRNA2<-unique(ceRNA_triples[,c(2,3)])
  colnames(mir_ceRNA2)<-c("mir","ceRNA2")
  ceRNA1_ceRNA2<-unique(ceRNA_triples[,c(1,3)])
  colnames(ceRNA1_ceRNA2)<-c("ceRNA1","ceRNA2")
  
##   
  mir_ceRNA1_c_p<-corrs(mir_ceRNA1,expr1=gene_exp,expr2=mir_M)
  mir_ceRNA2_c_p<-corrs(mir_ceRNA2,expr1=mir_M,expr2=gene_exp)
  ceRNA1_ceRNA2_c_p<-corrs(ceRNA1_ceRNA2,expr1=gene_exp,expr2=gene_exp)
  save(mir_ceRNA1_c_p,mir_ceRNA2_c_p,ceRNA1_ceRNA2_c_p,file="c_p.Rdata") 
 
#### 

   mir_ceRNA1_c_p_sig<-pair_filter(mir_ceRNA1_c_p,"mir","fdr",cthr1,pthr)
   rownames(mir_ceRNA1_c_p_sig)<-paste(mir_ceRNA1_c_p_sig[,1],mir_ceRNA1_c_p_sig[,2],sep="_")##以显著对命行名

   mir_ceRNA2_c_p_sig<-pair_filter(mir_ceRNA2_c_p,"mir","fdr",cthr2,pthr)
   rownames(mir_ceRNA2_c_p_sig)<-paste(mir_ceRNA2_c_p_sig[,1],mir_ceRNA2_c_p_sig[,2],sep="_")##以显著对命行名

   ceRNA1_ceRNA2_c_p_sig<-pair_filter(ceRNA1_ceRNA2_c_p,"gene","fdr",cthr3,pthr)
   rownames(ceRNA1_ceRNA2_c_p_sig)<-paste(ceRNA1_ceRNA2_c_p_sig[,1],ceRNA1_ceRNA2_c_p_sig[,2],sep="_")##以显著对命行名

   save(mir_ceRNA1_c_p_sig,mir_ceRNA2_c_p_sig,ceRNA1_ceRNA2_c_p_sig,file="sig_c_p.Rdata")
#############################################################################################
##### 

   miRNA_lncRNA_gene1<-filter1(ceRNA_triples,mir_ceRNA1_c_p_sig,"12")
   miRNA_lncRNA_gene2<-filter1(miRNA_lncRNA_gene1,mir_ceRNA2_c_p_sig,"23")
   sig_miRNA_lncRNA_gene<-filter1(miRNA_lncRNA_gene2,ceRNA1_ceRNA2_c_p_sig,"13")
   save(sig_miRNA_lncRNA_gene,file="sig_miRNA_lncRNA_gene.Rdata")
###############################################################################################
##
   
   cp1<-extract_cp(sig_miRNA_lncRNA_gene,"12",mir_ceRNA1_c_p_sig)
   cp2<-extract_cp(cp1,"13",ceRNA1_ceRNA2_c_p_sig)
   sig_miRNA_lncRNA_gene_cp<-extract_cp(cp2,"23",mir_ceRNA2_c_p_sig)
   save(sig_miRNA_lncRNA_gene_cp,file="sig_miRNA_lncRNA_gene_cp.Rdata")
   return(sig_miRNA_lncRNA_gene)
   
}

##
cor_alls<-function(a_flag,exprs,ceRNA_triples,ceRNA_mir_list,mir_ceRNA_list,cthr1,cthr2,cthr3,pthr){
     dir.create(a_flag) 
	 temp_path<-getwd();
	 temp_path<-paste(temp_path,"/",a_flag,sep="")
	 setwd(temp_path) 
	
	 a_list<-exprs[[a_flag]]
     gene_exp<-a_list[[1]]
	 
	 mir_M<-a_list[[2]]
	 result<-cor_all(ceRNA_triples,mir_M,gene_exp,cthr1,cthr2,cthr3,pthr)
	
    setwd("../") 
	return(result)
}

ceRNAs_sig<-function(geneid,cnv_M,mir_M,gene_exp, min_thr,ceRNA_mir_list,mir_ceRNA_list,cthr1,cthr2,cthr3,pthr){
## 
    dir.create(geneid,mode = "775") 
	temp_path<-getwd();
	temp_path<-paste(temp_path,"/",geneid,sep="")
	setwd(temp_path) 
    ceRNA_triples<-ceRNA_mir_ceRNA(geneid,ceRNA_mir_list,mir_ceRNA_list) 
	save(ceRNA_triples,file="ceRNA_triples.Rdata")
	flag_v<-cnv_M[geneid,]
	exprs<-case_normal_profile(flag_v,gene_exp,mir_M,min_thr)  
	cnv_flag<-names(exprs)
	
	sig_ceRNAs<-lapply(cnv_flag,cor_alls,exprs=exprs,ceRNA_triples=ceRNA_triples,ceRNA_mir_list=ceRNA_mir_list,mir_ceRNA_list=mir_ceRNA_list,cthr1=cthr1,cthr2=cthr2,cthr3=cthr3,pthr=pthr)
	 
	names(sig_ceRNAs)<-cnv_flag
	save(sig_ceRNAs,file="sig_ceRNAs.Rdata")
	setwd("../")
	return(sig_ceRNAs)


}
###################################################################
 
driver_lnc<-function(Dname,input,freq_thr1,exp_thr,freq_thr2,p_thr,min_thr,cthr1,cthr2,cthr3){
 


 
   path_save<-paste(input,"/",Dname,sep="")
   setwd(path_save)
   load("cnv_M.RData")
   load("lnc_M.Rdata")
   load("gene_M.Rdata")
   load("mir_M.Rdata")
   setwd(input)
   load("mir_ceRNA_list.Rdata")

 
   setwd(path_save)
   common_data(cnv_M,gene_M,lnc_M,mir_M)
   load("common_data.Rdata")
 
   cnv_M<-del_cnv(cnv_M, freq_thr1)
   
 
   lnc_M<-del_exp(lnc_M,exp_thr,freq_thr2)
 
   gene_M<-del_exp(gene_M,exp_thr,freq_thr2)
 
   mir_M<-del_exp(mir_M,exp_thr,freq_thr2)



   
   cnv_M<-common_data1(cnv_M,ceRNA_mir_list)
   gene_M<-common_data1(gene_M,ceRNA_mir_list)
   lnc_M<-common_data1(lnc_M,ceRNA_mir_list)
   mir_M<-common_data1(mir_M,mir_ceRNA_list)
   gene_exp<-rbind(lnc_M,gene_M)  
   
   mir_ceRNA_list<-mir_ceRNA_list[rownames(mir_M)] 
   ceRNA_mir_list<-ceRNA_mir_list[rownames(gene_exp)] 

   mir_ceRNA_list<-common_data2(gene_exp,mir_ceRNA_list) 
   ceRNA_mir_list<-common_data2(mir_M,ceRNA_mir_list)  
   
   cnv_M<-common_data1(cnv_M,ceRNA_mir_list) 
   gene_flag<-cnvs_exp(cnv_M,gene_exp,p_thr,min_thr) 
   cnv_M<-cnv_M[gene_flag,]
   
 
    save(cnv_M,file="cnv_M_for_driver.Rdata")
    genes<-rownames(cnv_M)
	result<-lapply(genes,ceRNAs_sig,cnv_M=cnv_M,mir_M=mir_M,gene_exp=gene_exp, min_thr=min_thr,ceRNA_mir_list=ceRNA_mir_list,mir_ceRNA_list=mir_ceRNA_list,cthr1=cthr1,cthr2=cthr2,cthr3=cthr3,pthr=p_thr)
    names(result)<-genes
	return(result)

}

#########################################################################################################################
#### 

is.null<-function(a_list){
       
	   if(length(unlist(a_list))==0){
	      return(T)
	   }else{
	      return(F)
	   }

}
extract_notnull<-function(all_sig_ceRNAs){
      T_F<-unlist(lapply(all_sig_ceRNAs,is.null))
	  result<-all_sig_ceRNAs[!T_F]
	  return(result)
}

##################################################################################
##  
remove_self1<-function(ceRNAs){
       if(length(ceRNAs)==0){
	      return(ceRNAs)
	   }else{
	      if(length(ceRNAs)==3){
		     return(NULL)
		  }
		  
	      pos<-which(ceRNAs[,1]==ceRNAs[,3])
		  if(length(pos)){
		     ceRNAs<-ceRNAs[-pos,]
		  }
		  return(ceRNAs)
	  }
}

remove_self2<-function(a_list){
    a_list<-lapply(a_list,remove_self1)
	return(a_list)
}
remove_self<-function(a_list){
    a_list<-lapply(a_list,remove_self2)
	return(a_list)
}
############################################################################################################

diff_net<-function(a_list){

    temp_l<-unlist(lapply(a_list,length))
	if(sum(temp_l==0)==0){ 
	   com_mir<-unique(intersect(a_list[[1]][,2],a_list[[2]][,2]))
	   if(length(com_mir)==0){ 
	     result<-rbind(a_list[[1]],a_list[[2]])
	   }else{ 
	      pos1<-which(a_list[[1]][,2]%in%a_list[[2]][,2])
		  pos2<-which(a_list[[1]][,3]%in%a_list[[2]][,3])
		  pos<-intersect(pos1,pos2)
		  if(length(pos)==0){ 
		    result<-rbind(a_list[[1]],a_list[[2]])
		  }else{ 
		    pos11<-which(a_list[[2]][,2]%in%a_list[[1]][,2])
		    pos22<-which(a_list[[2]][,3]%in%a_list[[1]][,3])
		    pos12<-intersect(pos1,pos2)
			
		    result<-rbind(a_list[[1]][-pos,],a_list[[2]][-pos12,])
		  }	   
	   }
	}else{ 
	  result<-a_list[temp_l!=0][[1]]
	}
	return(result)

}
