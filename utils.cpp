namespace Utils {
	void Swap(int &a, int &b){
		if (a == b) {
			return;
		}
		a ^= b;
		b ^= a;
		a ^= b;
		return;
	}
}

