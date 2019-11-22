#include "handle_database.h"
#include <cstring>

//Main prend un argument le chemin vers la database

void display_vector(const vector<char>* subject)
{
	for(unsigned int i=0; i<subject->size();++i)
	{
		if(i%80==0){cout << endl;}
		cout << subject->at(i) ;
	}
	
	cout << endl;
}

int compare_prot(int protlen,const string *protein,const vector<char>*dataprot){
	if (protlen ==1){
		if (protein->at(0)==dataprot->at(0))
			return 1;
		else
			return 0;
	}
	else{
		int nsize=(int) floor(protlen/2);
		string nprot1= protein->substr(0,nsize);
		string nprot2= protein->substr(nsize);
		vector<char> ndprot1;
		ndprot1.resize(nsize);
		for (int i=0; i<nsize; i++){
			ndprot1.at(i)=dataprot->at(i);
		}
		vector<char> ndprot2;
		ndprot2.resize(protlen-nsize);
		for (int i=nsize; i<protlen; i++){
			ndprot2.at(i-nsize)=dataprot->at(i);
		}
		if (compare_prot(nsize,&nprot1,&ndprot1) == 1 && compare_prot((protlen-nsize),&nprot2,&ndprot2) == 1)
			return 1;
		else
			return 0;
	}
}

int main(int argc, char* argv[])
{
	Handle_Database* database = new Handle_Database(argv[1]);
	
	string prot_sequence_query = "RRPARSGGDGGAPMTTGSRVVKYYDGSRRGSRRGLSTSGRSVKKDPAGLRDSLLSEDDRSAAAAPPPPPVHPVRDQLSSQLVRPSRGLGAYRTMSVFGSGWRPCRAAASHVRGAR";
	unsigned int length = prot_sequence_query.size();
	bool is_same = false ;
	bool is_different = false ;
	
	
	for(unsigned int i =0; i<database->get_database_size();++i)
	{
		if(database->get_size_sequence_prot(i) == prot_sequence_query.size())
		{
			/*
			database->fetch_prot_sequence(i); // on met a jour prot_active afin de contenir la i eme
			is_same = compare_prot(length, &prot_sequence_query, database->get_prot_active());
			*/
			is_different = false ;
			database->fetch_prot_sequence(i); // on met a jour prot_active afin de contenir la i eme
			for(unsigned int j=0; j<length; ++j)
			{
				if(database->get_prot_active()->at(j) != prot_sequence_query[j])
				{
					is_different = true ;
					break;
				}
			}
			
			if(!is_different)
			{
				cout << "Trouve index : " << i << endl ;
				break ;
			}
		}
	}
	
	delete database ;
	return 0 ;
}


