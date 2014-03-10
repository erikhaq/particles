#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

using namespace std;
int main( int argc, char **argv )
{
	vector<int> v;
	for(int i = 1; i < 10; i++)
	{
		int t = i*10;
		v.push_back(t);
		printf("T: %d\n", t);
	}

	std::vector<int>::iterator iter = v.begin();
	std::vector<int> v2;
	while(iter != v.end())
	{

		if(*iter == 50 || *iter ==30)
		{
			int t = *iter;
			v2.push_back(t);
			iter = v.erase(iter);
		}
		else{
			++iter;
		}
		
		printf("iter: %d\n", *iter);
		// ++iter;
	}
	iter = v2.begin();
	while (iter != v2.end())
	{
		printf("v2: %d\n", *iter);
		++iter;
	}
}