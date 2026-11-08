# Installed, offline SOFTWARE mock example; no genuine predictor or HLA proteins.
library(Neoantimon2)
root<-tempfile("indirect-class2-mock-");dir.create(root)
model<-file.path(root,"model");dir.create(file.path(model,"data"),recursive=TRUE)
f<-system.file("extdata","indirect-expanded","mock-class2.py",package="Neoantimon2")
stopifnot(file.copy(f,file.path(model,"netMHCIIpan")));Sys.chmod(file.path(model,"netMHCIIpan"),"0755")
writeLines("good",file.path(model,"data/mode.txt"))
m<-jsonlite::read_json(system.file("extdata","indirect-autorun","synthetic-reference.json",package="Neoantimon2"));template<-m$records[[1]]
alleles<-c("A*02:01","A*01:01","DRA*01:01","DRB1*01:01","DRB3*01:01","DRB4*01:01","DRB5*01:01")
m$records<-lapply(alleles,function(a){r<-template;r$allele<-paste0("HLA-",a);r$sequence_identifier<-paste0("SYNTHETIC:",a)
 r$protein_sequence<-if(a=="A*02:01")"ACDEFGHIKLMNPQRSTVWY"else"ACDEFGHIKAAAAAA";n<-nchar(r$protein_sequence)
 r$domains<-r$covered_regions<-list(list(name="full_precursor",start=1L,end=n));r})
ref<-file.path(root,"reference.json");jsonlite::write_json(m,ref,auto_unbox=TRUE,null="null",digits=NA)
run<-neo_indirect_class2(donor=list(A=rep("A*02:01",2)),
 recipient=list(A=rep("A*01:01",2),DRA=rep("DRA*01:01",2),DRB1="DRB1*01:01",DRB3="DRB3*01:01",DRB4="DRB4*01:01",DRB5="DRB5*01:01"),
 installation=list(path=model,authorization=list(status="mock_only",terms="Synthetic software demonstration; no genuine model")),
 reference=list(path=ref,sha256=digest::digest(file=ref,algo="sha256",serialize=FALSE)),
 source_loci="A",self_loci="A",output_dir=file.path(root,"run"))
stopifnot(run$status=="COMPLETED",!run$manifest$real_predictor_executed,run$result$group_summary[[1]]$complete_count==20L)
cat("MOCK SOFTWARE EXAMPLE ONLY; 20 artificial units; no biological or clinical validation.\n")
