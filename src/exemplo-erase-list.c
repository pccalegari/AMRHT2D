// erasing from list
#include <iostream>
#include <list>
using namespace std;

int main ()
{
  std::list<int> mylist;
  std::list<int>::iterator it1,it2;
  int vc, ct = 0;
  // set some values:
  for (int i=1; i<10; ++i) mylist.push_back(i*10);

                              // 10 20 30 40 50 60 70 80 90
  it1 = it2 = mylist.begin(); // ^^
  advance (it2,6);            // ^                 ^
  ++it1;                      //    ^              ^

  it1 = mylist.erase (it1);   // 10 30 40 50 60 70 80 90
                              //    ^           ^

  it2 = mylist.erase (it2);   // 10 30 40 50 60 80 90
                              //    ^           ^

  ++it1;                      //       ^        ^
  --it2;                       //       ^     ^
  
  mylist.insert(it1, 10);
  advance(it1,-2);
  it1 = mylist.erase(it1);
  
  //mylist.erase (it1,it2);     // 10 30 60 80 90
                              //        ^

  std::cout << "mylist contains:";
  for (it1=mylist.begin(); it1!=mylist.end(); ++it1)
    std::cout << ' ' << *it1;
  std::cout << '\n';
  cin >> vc;
  for (it1=mylist.begin(); it1!=mylist.end(); ++it1){
    //cout << vc << " " << *it1 << endl;
    if(vc < *it1){
      mylist.insert(it1,vc);
      cout << "A" << endl;
      break;
    }
    ct++;
  }
  cout << ct << " " << mylist.size() << endl;
  if(ct == mylist.size())
    mylist.insert(mylist.end(),vc);
  
  for (it1=mylist.begin(); it1!=mylist.end(); ++it1)
    std::cout << ' ' << *it1;
  
  std::cout << '\n';
  
  return 0;
}
