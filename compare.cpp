int compare_prot(int protlen, string *protein, vector<char>*dataprot){
	if (protlen ==1){
		if (*prot==*dataprot[0])
			return 1;
		else
			return 0;
	}
	else{
		int nsize=(int) floor(protlen/2);
		string nprot1= prot->substr(0,nsize-1);
		string nprot2= prot->substr(nsize);
		vector<char> ndprot1;
		for (int i=0, i<nsize, i++)
			ndprot1[i]=(*dataprot)[i];
		vector<char> ndprot2; 
		for (int i=nsize; i<protlen; i++)
			ndprot2[i]=(*dataprot)[i-nsize];
		if (compare_prot(nsize,&nprot1,&ndprot1) == 1 && compare_prot(nsize+1,&nprot2,&ndprot2) == 1)
			return 1;
		else
			return 0;
	}
}
