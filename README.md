# AnalysisiRPC


## Github staff
1 ) Fork the repository  
2 ) git remote add <your repo name> <your repo link.git>  
3 ) git fetch <your repo name>  
4 ) git checkout -b localmaster <your repo name>/main  
After this stage you should have all the scripts that are available on the gihub.  
When you change things, please firt commit as below:  
5 ) git commit -am "Write describing comments"  
6 ) git push <your repo name> localmaster  

Then in your we browser do a pull request firt to your main branch accept it merge it.   
Check your changes, if everything is OK,then make a pull request to Ece's repository. 
 
## Analysis

### Where is the data ?

We have two data types: raw (.dat) and ROOT (.root) version  
Depending on the data taking location our raw and root files are stored in different folders.  
At the moment we are taking data in Gif++, we are using the machine in the preparation area.  
You can access the machine (acqcmsmu01) via ssh: ssh acqcms@128.141.151.164  

Raw files are stored in this direcotry:/ 
/data/gifOctober/raw/  

1 ) If you don't have the root fie for necessary run, please produce it.  
Production:  

