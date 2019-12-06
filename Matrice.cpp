map<char,int> blos_dictionnary;
blos_dictionnary= { {'A',0}, {'R',1}, {'N',2}, {'D',3}, {'C',4}, {'Q',5}, {'E',6}, {'G',7}, {'H',8}, {'I',9}, {'L',10}, {'K',11}, {'M',12}, {'F',13}, {'P',14}, {'S',15}, {'T',16}, {'W',17}, {'Y',18}, {'V',19}, {'B',20}, {'Z',21}, {'X',22},{'*',23}};

int matrice(string Res1, string Res2, const map<char,int>* blos_dictionnary, const vector<vector<int>>* matri){
	i =0;
	res=0;

	while(Res1[i] != ' ' && Res2[i]!= ' '){
		res+=matri->at((*blos_dictionnary)[Res1[i]])[(*blos_dictionnary)[Res2[i]]);
		i++;
	}
	cout <<"Nombre de caractères: " << i << " Résultat: " << res << endl;
	return res;
}
