#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <map>

using namespace std;

#include "constants.h"


#include "grammar_reader.h"

//********************** The Code for Grammer Reading and Validation Begins ********
unsigned long SCFGrule::NTcount = 0;
unsigned long SCFGrule::Tcount = 0;
unsigned long SCFGrule::RuleCount = 0;
string SCFGrule::NT;
string SCFGrule::T;
string SCFGrule::NTonLHS;

//Function to initialize rule
void SCFGrule::InitializeRule(double & score, string & left, string & right)
{
  RuleScore     = score;
  ruleLHS       = left;
  ruleRHS       = right;
  RuleNumber        = RuleCount+1;
  RuleType      = _UNKNOWN_RULE;
  UsedOnceFlag      = false;

  No_of_RHSsymbols = 0;
}

//Function to set OnceUsedFlag if rule has been used at least once
void SCFGrule::SetUsedOnceFlag(void)
{
    UsedOnceFlag = true;
}

// Function for rule checking the grammar and update parameters has to be run for each object, optional
//optional
void SCFGrule::CountUpdate(void)
{
    string dummy, newRHS;
    string star = "*";
    dummy = star + ruleLHS + star;
    //Updating each terminal or non-terminal by enclosing it in *X*
    

    if(!isupper(ruleLHS[0]))
    {
        cout << endl << "Lowercase Symbol on LHS : " << ruleLHS 
             << " at Rule No. " << RuleNumber << endl;
        exit(1);
    }

    ruleLHS = dummy;

    if(NT.empty())
    {   
        NT = dummy; 
        ++NTcount; 
        NTonLHS = NT;
    } 
    else if(NT.find(dummy) == string::npos) //dummy not found
    { 
        NT += dummy; 
        NTonLHS += dummy; 
        ++NTcount; 
    }
    else if(NTonLHS.find(dummy) == string::npos) 
    {
        NTonLHS += dummy; 
    }

    istringstream iss(ruleRHS);
    string temp;
    
    int count = 0;
    //This variable keeps track of how many symbols read from RHS
     

    while (iss >> temp)
    {
        ++count;
        if(!isupper(temp[0]))
        {//If the rule does not begin with upper case
            if((count != 1) || (!isalpha(temp[0])))
            {
                cout << endl << "Lowercase Symbol on RHS : " << ruleRHS 
                     << " at Rule No. " << RuleNumber << endl;
                exit(1);
            }
            else 
            {
                RuleType = _TERM_RULE;
                dummy = star + temp + star;
                if(T.empty())
                {   T = dummy; ++Tcount; } 
                else if((T.find(dummy) < 0) || (T.find(dummy) >= T.length()))
                { T += dummy; ++Tcount; }

                newRHS += dummy;

            }
            break;

        }

        else
        {//If the rule begins with upper case= it is a non-terminal rule
            dummy = star + temp + star;
        
            //int x = NT.find(dummy);
            RuleType = _NONTERM_RULE;

            if(NT.find(dummy) == string::npos)
            {NT += dummy;
                        ++NTcount;//To confirm from Vipin- Was the updation of NTcount deliberately avoided here?
                        }
            newRHS += dummy;
        }
        temp.erase(temp.begin(),temp.end());
    }
    ruleRHS = newRHS;

    No_of_RHSsymbols = count;
    ++RuleCount;
}

// validation of Grammar has to be run once (optional)

bool SCFGrule::ValidateGrammar(void)
{
    //validate Terminal strings
    for(int i = 0; i < T.length() ; ++i)
    {
        if((T[i] == '*') || (islower(T[i]))) continue;
        else
        {
            cout << endl <<"Some Error in Terminal string : " << T.substr(i-1,3) << endl << T;
            return false;
        }
    }                                 

    //validation of Non Terminal String
    //Checking if all NTs occur on LHS
/*
    if(NT.length() != NTonLHS.length())
    {
        cout << endl << "NT and NTonLHS are different in size" << endl
             << "NT      : " << NT << endl
             << "NTonLHS : " << NTonLHS << endl;
        return false;
    }
*/
    //char* temp = new char[NT.length() + 1];
    
    string sortedNT = NT;
    string sortedNTonLHS = NTonLHS;

    sort(sortedNT.begin(),sortedNT.end());
    sort(sortedNTonLHS.begin(),sortedNTonLHS.end());

    if(sortedNT != sortedNTonLHS)
    {
        cout << endl << "NT and NTonLHS are different in nature" << endl
             << "NT      : " << NT << endl
             << "NTonLHS : " << NTonLHS << endl
             << "T       : " << T << endl;
        return false;
    }

    return true;
}


// This function reads the grammar and assumes the grammar according to following format
// Each Rule on a separate line : <probability> <Rule LHS> --> <Rule RHS>
// Each LHS can have a single NT name with no space as space is a delimiter for NTs
// RHS should have either a T or string of NTs separated by a space
// Also small letters are assumed to be Terminals (T) and Capital ones as Non-Terminals (NT)

vector<SCFGrule> ReadGrammarFile(ifstream * grammarfileptr)
{
vector<SCFGrule> TempRuleList;

    if(grammarfileptr == NULL)
    {
        grammarfileptr = new ifstream;
        if(grammarfileptr == NULL)
        {
            cout << endl << "Unable to Read or Allocate memory for reading Grammar File" << endl;
            exit(1);
        }
        char filename[MAXPATH+1];
        cout << endl << "\n Enter Name of Grammar File : ";
        cin.getline(filename,MAXPATH);
            grammarfileptr->open(filename, ios::in);
        if(!grammarfileptr->good())
        {
            cout << endl << "Unable to Read Grammar File : " << filename << endl;
            exit(1);
        }
    }

    // Start Parsing the Grammar file

        SCFGrule tempruleobj;

    while(!grammarfileptr->eof())
    {
        char temprule[MAXPATH];
        grammarfileptr->getline(temprule,MAXPATH-1);
        if(!strlen(temprule)) continue;
        
        double score;
        istringstream iss(temprule);
        string rule;

        if((iss >> score) && (score >= 0.0) && (score <= 1.0));
        else
        {
            cout << "\n Error Reading score in rule !!\n"
                 << temprule;
            exit(1);
        }

        rule = iss.str();
        int tempitr = -1;
        while(!isspace(rule[++tempitr]));
        rule.erase(0,tempitr);
                tempitr = -1;
        while(!isalpha(rule[++tempitr]));
        rule.erase(0,tempitr);
        //this is done because to remove the initial prob value from input line
        int delimiter = rule.find("-->");

        if((delimiter < 0) || (delimiter >= rule.length()))
        {
            cout << "\n Error Reading --> in rule !!\n"
                 << temprule;
            exit(1);
        }
        
        string tempLHS, tempRHS;
        tempLHS = rule.substr(0,delimiter);
        tempRHS = rule.substr(delimiter+3,rule.length()-delimiter-3);

        //PRUNING STARTING AND TRAILING WHITESPACES
        long spacecount = -1;
        while(!isalpha(tempLHS[++spacecount]));
        tempLHS.erase(0,spacecount);

        spacecount = tempLHS.length();
        while(!isalnum(tempLHS[--spacecount])); 
        tempLHS.erase(spacecount+1, tempLHS.length() - spacecount -1);

        spacecount = -1;
        while(!isalpha(tempRHS[++spacecount]));
        tempRHS.erase(0,spacecount);

        spacecount = tempRHS.length();
        while(!isalnum(tempRHS[--spacecount])); 
        tempRHS.erase(spacecount+1, tempRHS.length() - spacecount -1);
        //PRUNING ENDS

        tempruleobj.InitializeRule(score, tempLHS, tempRHS);
        tempruleobj.CountUpdate();
                
        TempRuleList.push_back(tempruleobj);
    }

    //IF OPTIONAL ERRORCHECK REQUIRED PUT HERE
    if(!SCFGrule::ValidateGrammar())
    {
        cout << "\n Error in Grammar.\n";
        exit(1);
    }
return (TempRuleList);
}
//********************** The Code for Grammer Reading and Validation Ends **********
