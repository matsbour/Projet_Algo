#include "myProtein.h" 


myProtein::myProtein(const string proteinfile)
{
	header = output_header(proteinfile);
	sequence = output_sequence(proteinfile);
	size=sequence.size();
}
	

// Getters
	 
string* myProtein::getSequence(){return &sequence;}
string* myProtein::getHeader(){return &header;}
const int myProtein::getSize() const{return size;}
	 

string myProtein::output_header(const string proteinfile) 
{ 
	/**
	 * @desc extrait la description de la protéine 
	 * @param string : fichier FASTA de la protéine
	 * @return name : première ligne du fichier **/
   
   ifstream file(proteinfile,std::ifstream::binary);
	if(file.is_open()){
    
		string line, name;
		while(getline(file, line)) 
			{ 
			
			if(line[0] == '>'){ // marqueur identifiant la première ligne
			  name=line;  
			  break;
			} 
		}
	  return name;
	}
	
	else{cout << "Cannot read: " << proteinfile<<endl;} 
	file.close();
	
	return ""; //Retour vide = error
}
	

string myProtein::output_sequence(const string proteinfile) 
{ 
   /**
	 * @desc extrait la séquence de la protéine 
	 * @param string : fichier FASTA de la protéine
	 * @return name : séquence de la protéine **/
   
    ifstream file(proteinfile,std::ifstream::binary);
	if(file.is_open()){

    
    string content,line;
	while(getline(file, line)) 
		{ 
		if(line[0] != '>'){ // vérifie que c'est pas la première ligne
			content=content+line;} 
	}
	
  return content;}
  
  else{cout << "Cannot read: " << proteinfile<<endl;}
	file.close();
	
	return "";//Vide = error
}

