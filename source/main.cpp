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

int main(int argc, char* argv[])
{
	Handle_Database* database = new Handle_Database(argv[1]);
	
	string prot_sequence_query = "RRPARSGGDGGAPMTTGSRVVKYYDGSRRGSRRGLSTSGRSVKKDPAGLRDSLLSEDDRSAAAAPPPPPVHPVRDQLSSQLVRPSRGLGAYRTMSVFGSGWRPCRAAASHVRGAR";
	unsigned int length = prot_sequence_query.size();
	bool is_different = false ;
	
	
	for(unsigned int i =0; i<database->get_database_size();++i)
	{
		if(database->get_size_sequence_prot(i) == prot_sequence_query.size())
		{
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


