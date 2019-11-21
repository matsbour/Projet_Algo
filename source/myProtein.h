#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class myProtein {

	public:
	
	myProtein();
	
	// Instance variables 
	string header;
	string sequence;
	
	// Function declaration of output_header() to extract header from file Data Base 
	string output_header();
	
	// Function declaration of output_sequence() to extract sequence from file Data Base 
	string output_sequence();

};
