#include "handle_query.h"

Handle_query::Handle_query(string filepath, string database_filepath)
{
	this->database_filepath = database_filepath ;
	
	ifstream finput (filepath) ;
	
	if(finput.is_open())
	{
		this->query_description = finput.getLine() ;
		//Traiter la string query pour enlever > du debut
		string temp = "" ;
		while( (temp=finput.getLine()) != EOF ) // tant que nous n avons pas atteint la fin on lit ligne par ligne
		{
			this->query_prot_chain+=temp ;//On forme la liste complete de protein en concatenant les temp
		}
	}
	else
	{
		exit 1 ;
	}
	
	finput.close() ;
}

string get_prot_chain() const{return this->query_prot_chain;}
string get_query_description() const{return this->query_description;}

