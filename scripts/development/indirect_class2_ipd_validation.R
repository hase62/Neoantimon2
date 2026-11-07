# Opt-in genuine technical E2E. No patient data; no bundled IPD/model assets.
# Rscript --vanilla scripts/development/indirect_class2_ipd_validation.R INSTALL CACHE NEW_OUTPUT
args<-commandArgs(TRUE)
stopifnot(length(args)==3L)
library(Neoantimon2)
install<-normalizePath(args[[1]],mustWork=TRUE);cache<-args[[2]];dest<-args[[3]]
stopifnot(!dir.exists(dest));dir.create(dest,recursive=TRUE)
auth<-list(status="authorized_local_installation",terms="Explicit local execution of user-acquired licensed models for academic technical validation; no redistribution")
base<-list(installation=list(path=install,authorization=auth),reference=list(release="3.66.0",cache=cache),batch_size=1024L)
absent<-list(status="absent",copy_number=0)
rec<-list(A=rep("A*01:01",2),DRA=rep("DRA*01:01",2),DRB1=rep("DRB1*01:01",2),DRB3=absent,DRB4=absent,DRB5=absent)
pair<-function(b)list(list(pair=paste0("DRA*01:01~",b),pairing="explicit_expressed",expression_confirmed=TRUE,evidence="Explicit hypothetical research declaration; not observed patient expression"))
common<-list(A=rep("A*02:01",2),B=rep("B*07:02",2),C=rep("C*07:02",2),DRB1=rep("DRB1*03:01",2),DRB3="DRB3*01:01",DRB4=absent,DRB5=absent,
 DQA1=rep("DQA1*01:01",2),DQB1=rep("DQB1*05:01",2),DPA1=rep("DPA1*02:02",2),DPB1=rep("DPB1*19:01",2),MICA=rep("MICA*001",2),MICB=rep("MICB*002",2))
recipient_common<-common;recipient_common$A<-rep("A*01:01",2);recipient_common$DRA<-rep("DRA*01:01",2)
cases<-list(DRB1=list(donor=list(A=rep("A*02:01",2)),recipient=rec,profile="dr-expanded-v1",source_loci="A",self_loci="A"),
 DRB3_all13=list(donor=common,recipient=recipient_common,profile="dr-expanded-v1",presenters=pair("DRB3*01:01")),
 DRB4=list(donor=list(A=rep("A*02:01",2),DRB4="DRB4*01:01"),recipient=list(A=rep("A*01:01",2),DRA=rep("DRA*01:01",2),DRB4="DRB4*01:01"),profile="dr-expanded-v1",presenters=pair("DRB4*01:01"),source_loci=c("A","DRB4"),self_loci=c("A","DRB4")),
 DRB5=list(donor=list(A=rep("A*02:01",2),DRB5="DRB5*01:01"),recipient=list(A=rep("A*01:01",2),DRA=rep("DRA*01:01",2),DRB5="DRB5*01:01"),profile="dr-expanded-v1",presenters=pair("DRB5*01:01"),source_loci=c("A","DRB5"),self_loci=c("A","DRB5")),
 DQ_DP=list(donor=list(A=rep("A*02:01",2)),recipient=list(A=rep("A*01:01",2),DQA1=rep("DQA1*01:01",2),DQB1=rep("DQB1*05:01",2),DPA1=rep("DPA1*02:02",2),DPB1=rep("DPB1*19:01",2)),profile="dqdp-exploratory-v1",source_loci="A",self_loci="A",presenters=list(list(pair="DQA1*01:01~DQB1*05:01",pairing="confirmed_cis",expression_confirmed=TRUE,evidence="Fictional research pairing; not measured"),list(pair="DPA1*02:02~DPB1*19:01",pairing="explicit_expressed",expression_confirmed=TRUE,evidence="Fictional research pairing; not measured"))))
summary<-list()
for(name in names(cases)){
 a<-c(cases[[name]],base,list(output_dir=file.path(dest,name)))
 saveRDS(a,file.path(dest,paste0(name,"-args.rds")))
 cat("BEGIN",name,"\n");flush.console()
 x<-do.call(neo_indirect_class2,a)
 stopifnot(x$status=="COMPLETED",x$manifest$real_predictor_executed,x$manifest$provider_version=="4.3k")
 saveRDS(x,file.path(dest,paste0(name,"-result.rds")))
 summary[[name]]<-list(status=x$status,origins=x$result$donor_origin_count,self_origins=x$result$self_origin_count,queries=x$result$unique_query_count,
  candidates=length(x$result$candidates),links=length(x$result$candidate_origin_links),groups=x$result$group_summary)
 cat("COMPLETED",name,x$result$donor_origin_count,"origins\n");flush.console()
}
jsonlite::write_json(summary,file.path(dest,"summary.json"),auto_unbox=TRUE,null="null",pretty=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(dest,"sessionInfo.txt"))
cat("ALL_REAL_CASES_COMPLETED\n")
