#include "myProtein.h" 


//rajouter taille chaine de caractères

myProtein::myProtein()
{
	//header = " " ;
	header = output_header();
	sequence = output_sequence();
	size=sequence.size();
}
	

	 // Getters
	 
	
	string myProtein::getSequence() const{return sequence;}
	
	string myProtein::getHeader() const{return header;}

	int myProtein::getSize() const{return size;}
	 

// Function definition of output_header() to extract header from file Data Base 

string myProtein::output_header() 
{ 
   // Object to read from file 
    ifstream in("proteine.fasta");
    
    string line, name;
	while(getline(in, line)) 
		{ 
		
		if(line[0] == '>'){ // Identifier marker
          name=line;  
          
          //std::cout << "name:" <<name << std::endl; 
                
          break;
		} 
	}
	
	in.close();
  
  return name;
} 

// Function definition of output_sequence() to extract sequence from file Data Base 

string myProtein::output_sequence() 
{ 
   // Object to read from file 
    ifstream in("proteine.fasta"); 

    
    string content,line;
	while(getline(in, line)) 
		{ 
		
		if(line[0] != '>'){ // Identifier marker
			content=content+line;  
				
			//std::cout << " content:" <<content  << std::endl;    
          
		} 
	}
	
  return content;
}
