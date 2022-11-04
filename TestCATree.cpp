#include <iostream>
#include <vector>
using namespace std;

#include "CATree.h"

int main()
{
   // Declare the node object vectors
   vector<Node *> NodesEScheme, NodesWTAScheme, NodesCATree;

   // Loop over jet constituents.  Here we use a toy with 25 particles
   for(int i = 0; i < 25; i++)
   {
      // Set particle kinematics
      FourVector P;
      P.SetPtEtaPhiMass(10 + i * 0.0001, (i / 5) * 0.01, (i % 5) * 0.01, 0);

      // Add into the node object vector
      NodesEScheme.push_back(new Node(P));
      NodesWTAScheme.push_back(new Node(P));
      NodesCATree.push_back(new Node(P));
   }

   // Do the reclustering!
   BuildCATree(NodesEScheme, -1, EScheme);
   BuildCATree(NodesWTAScheme, -1, WTAScheme);
   BuildCATree(NodesCATree, 0, EScheme);

   // Some output
   cout << "EScheme: " << NodesEScheme[0]->P.GetEta() << " " << NodesEScheme[0]->P.GetPhi() << endl;
   cout << "WTAScheme: " << NodesWTAScheme[0]->P.GetEta() << " " << NodesWTAScheme[0]->P.GetPhi() << endl;

   // Find soft drop node.  In this case z = 0.1, beta = 0.0, RJet = 0.4
   Node *SDNode = FindSDNode(NodesCATree[0], 0.1, 0.0, 0.4);
   if(SDNode->N > 1)
   {
      FourVector Subjet1 = SDNode->Child1->P;
      FourVector Subjet2 = SDNode->Child2->P;
      double ZG = min(Subjet1.GetPT(), Subjet2.GetPT()) / (Subjet1.GetPT() + Subjet2.GetPT());
      double RG = GetDR(Subjet1, Subjet2);
      double MG = (Subjet1 + Subjet2).GetMass();
      cout << "Groomed result (0.1, 0.0): " << ZG << " " << RG << " " << MG << endl;
   }
   else
      cout << "Groomed result (0.1, 0.0): -1 -1 -1" << endl;
   
   // Do another soft drop setting
   SDNode = FindSDNode(NodesCATree[0], 0.5, 1.5, 0.4);
   if(SDNode->N > 1)
   {
      FourVector Subjet1 = SDNode->Child1->P;
      FourVector Subjet2 = SDNode->Child2->P;
      double ZG = min(Subjet1.GetPT(), Subjet2.GetPT()) / (Subjet1.GetPT() + Subjet2.GetPT());
      double RG = GetDR(Subjet1, Subjet2);
      double MG = (Subjet1 + Subjet2).GetMass();
      cout << "Groomed result (0.5, 1.5): " << ZG << " " << RG << " " << MG << endl;
   }
   else cout << "Groomed result (0.5, 1.5): -1 -1 -1" << endl;

   // Clean up is important.  We only need to delete the root node and everything will be cleaned
   delete NodesEScheme[0], NodesWTAScheme[0], NodesCATree[0];

   return 0;
}





