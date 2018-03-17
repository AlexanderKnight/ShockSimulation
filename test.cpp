#include <iostream>
using namespace std;
#include <vector>

int main ()
{
	for (int i = 0; i <=10; i++){
		cout << i << endl;
	}
	std::vector<int> test (10,3);

	for (int i = 0; i <= 15; i++){
		cout << "i=" << i << ", test[i]=" << test[i] << endl;
	}
	cout << test[10] << endl;
	return 0;
}

