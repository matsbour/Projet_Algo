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
	 



// Function definition of output_header() to extract header from file Data Base 

string myProtein::output_header(const string proteinfile) 
{ 
   // Object to read from file 
   
   ifstream file(proteinfile,std::ifstream::binary);
	if(file.is_open()){
    
		string line, name;
		while(getline(file, line)) 
			{ 
			
			if(line[0] == '>'){ // Identifier marker
			  name=line;  
			  break;
			} 
		}
		
	  return name;
	}
	
	else{cout << "Cannot read: " << proteinfile<<endl;}
	file.close();
	
	return ""; //Retou vide = error
}
	
 

// Function definition of output_sequence() to extract sequence from file Data Base 

string myProtein::output_sequence(const string proteinfile) 
{ 
   // Object to read from file 
   
    ifstream file(proteinfile,std::ifstream::binary);
	if(file.is_open()){

    
    string content,line;
	while(getline(file, line)) 
		{ 
		if(line[0] != '>'){ // Identifier marker
			content=content+line;} 
	}
	
  return content;}
  
  else{cout << "Cannot read: " << proteinfile<<endl;}
	file.close();
	
	return "";//Vide error
}

