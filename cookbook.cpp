//math-1
	//modular airthmatic
	//**(a+b)%m : (a%m + b%m)%m;
	//**(a-b)%m : (a%m - b%m)%m;
	//**(a*b)%m : (a%m * b%m)%m;
	//**(a/b)%m : (a%m - b'%m)%m; b' is modular inverse of b under m
	//exponential power O(log(y)) ** we are calculating x^y
	long long power(long long x,long long y,long long mod)
	{
		long long res=1;
		while(y){
			if(y%2==1)
			res=(res*x)%mod;
			x=(x*x)%mod;
			y/=2;	
		}
		return res;
	}
	//find gcd gcd(a,b)=gcd(b,a%b) [euclid theorm] O(log(max(a,b)))
	long long gcd(long long a,long long b)
	{
		if(b==0)
			return a;
		else
			return gcd(b,a%b);
	}
	//extended euclid thoerm : for both gcd and multiplicative inverse O(log(max(a,b)))
	long long d, x, y;
	void extendedEuclid(long long  A,long long B) {
    	if(B == 0) {
        	d = A;
        	x = 1;
        	y = 0;
    	}
    	else {
    	    extendedEuclid(B, A%B);
    	    int temp = x;
    	    x = y;
    	    y = temp - (A/B)*y;
    	}
	}
	//modular inverse using extended euclid theorm
	long long modinverse(long long a,long long m)
	{
		extendedEuclid(a,m);
		x=(x+m)%m;
		return x;
	}
	//modular inverse using fermat's little theorm **can be only implied when m is prime**
	long long modinverse1(long long a,long long m)
	{
		return power(a,m-2,m);
	}
//math-2
	//check whether prime or not O(sqrt(n))
	bool isprime(long long n)
	{	
		if(n==1)
			return false;
		for(long long i=2;i*i<=n;i++)
			if(n%i==0)
				return false;
		return true;
	}
	//find all prime number upto a number N. O(NloglogN)[sieve of eranthoses] N~10^6
	isprime[n+1]={0}//denotes every number is prime
	void sieve(long long n)
	{
		isprime[0]=isprime[1]=1;
		for(int i=2;i*i<=n;i++)
			if(isprime[i]==0)
			for(int j=2*i;j<=n;j+=i)// j can start from i*i as number less than i must have already marked 'j' as non prime
				isprime[j]=1;
		//number which are prime are now marked with 0 in array while not prime are marked by 1
			for(int i=2;i<=n;i++)
				if(isprime[i]==0)
					cout<<i<<' ';
	}
	//factorize numbers O(sqrt(n))
	vector<int> factorize(int n)
	{	vector<int> res;
		for(int i=2;i*i<=n;i++)
		{
			while(n%i==0)
			{
				res.push_back(i);
				n/=i;
			}
		}
		if(n!=1)
			res.push_back(n);
		return res;
	}
	//factorization can be made fast using sieve to keep track of minimum prime factor of a number
	int minPrime[n + 1];
	for (int i = 2; i * i <= n; ++i) {
	    if (minPrime[i] == 0) {         //If i is prime
	        for (int j = i * i; j <= n; j += i) {
	            if (minPrime[j] == 0) {
	                minPrime[j] = i;
	            }
	        }
	    }
	}
	for (int i = 2; i <= n; ++i) {
	    if (minPrime[i] == 0) {
	        minPrime[i] = i;
	    }
	}
	vector<int> factorize(int n) {
    vector<int> res;
    while (n != 1) {
        res.push_back(minPrime[n]);
        n /= minPrime[n];
    }
    return res;
	}
//searching algorithm(binary-search) O(logN) 
	//to find a number just smaller than a number in an array
	long lower_bound(long *a,long p,long start,long end)
	{   
		if(p<a[0])
    	return -1;
    	long mid;
    	while(start<end)
    	{
        	mid=ceil((start+end)/2.0);
        	if(a[mid]==p)
        	return a[mid];
        	else if(a[mid]>p)
        	{
        	    end=mid-1;
        	}
        	else{
        	    start=mid;
        	}
        	mid=ceil((start+end)/2.0);
    	}
    	return a[mid];
	}
	//to find a number just greater than a given number from an array
	long upper_bound(long *a,long p,long start,long end)
	{   
		if(p>a[end])
	    	return -1;
	    long mid;
	    while(start<end)
	    	{
	    	    mid=(start+end)/2;
	    	    if(a[mid]<p)
	    	    {
	    	        start=mid+1;
	    	    }else end=mid;
	    	    mid=(start+end)/2;
	    	}
	    return a[mid];
	}	
//Bit-Manipulation
	//count no_of_one O(no_of_one_in_binary_equivalent_of_number) : worst case- O(logN)
	int count_one(long n)
	{
		int count=1;
		while(n)
		{
			n=n&(n-1);
			count++;
		}
		return count;
	}
	//check if i-th bit is set or not
	bool check (int N)
    {
        if( N & (1 << i) )
            return true;
        else
            return false;
    }
//Dynamic-Programming
    //*******************************************
    //DP - Length of Longest common subsequence
    //*******************************************
    string first,second;
    int LCS(int lf,int rf,int ls,int rs)
    {
    	if(lf>rf||ls>rs)
    		return 0;
    	if(first[rf]==second[rs])
    		return 1+LCS(lf,rf-1,ls,rs-1);
    	else return max(LCS(lf,rf,ls,rs-1),LCS(lf,rf-1,ls,rs));
    }
    void LCSDriver()
    {
    	cin>>first>>second;
    	cout<<LCS(0,first.length()-1,0,second.length()-1);
    }
    //*************************************
    //DP - Longest Increasing subsequence
    //*************************************
    int LIS(int length,int *a)
    {
    	for(int i=0;i<length;i++)
    		lis[i]=1;
    	for(int i=1;i<length;i++)
    		for(int j=0;j<i;j++)
    		{
    			if(a[i]>a[j])
    			{
    				if(lis[i]<=lis[j])
    					lis[i]=lis[j]+1;
    			}
    		}
    	int ans=-1;
    	for(int i=0;i<length;i++)
    		if(ans<lis[i])
    			ans=lis[i];
    	return ans;
    }
    //***************************************************************************************
    //Edit Distance - given two strings find out minimum cost to convert string 1 to string 2 
    //possible operation -insert,replace,remove
    //len of 1st string=m : len of 2nd string=n
    //****************************************************************************************
    int min(int a,int b,int c)
    {
    	return min(min(a,b),c);
    }
    int editDistance(string st1,string st2,int m,int n)
    {
    	int dp[m+1][n+1];
    	for(int i=0;i<=m;i++)
    	{
    		for(int j=0;j<=n;j++)
    		{
    			if(i==0)
    				dp[i][j]=j;
    			else if(j==0)
    				dp[i][j]=i;
    			else{
    				if(st1[i-1]==st2[j-1])
    					dp[i][j]=dp[i-1][j-1];
    				else dp[i][j]=1+min(dp[i-1][j-1],dp[i-1][j],dp[i][j-1]);
    			}
    		}
    	}
    	return dp[m][n];
    }
    //************************
    //DP-0-1 knapsack problem
    //************************
    int knapsack(int n,int e)
	{
		int val[n+1],cost[n+1];
		for(int i=1;i<=n;i++)
			cin>>val[i];
		for(int i=1;i<=n;i++)
			cin>>cost[i];
		int dp[e+1][n+1];
		memset(dp,0,sizeof(dp));
		for(int i=1;i<=e;i++)
		{
			for(int j=1;j<=n;j++)
			{
				if(cost[j]<=i)
					dp[i][j]=max((val[j]+dp[i-cost[j]][j-1]),dp[i][j-1]);
				else dp[i][j]=dp[i][j-1];
			}
		}
	//	cout<<dp[e][n];
		return dp[e][n];
	}
	//****************************************************************
	//DP-partition problem
	// i.e. find if array can be divided into two subset of equal sum
	//****************************************************************
	bool partition(int *a, int n,int sum)
	{
		if(sum%2)
			return false;

		int dp[sum/2 +1 ][n+1];
		for(int i=0;i<n+1;i++)
			dp[0][i]=true;
		for(int i=1;i<=sum/2 ;i++)
			dp[i][0]=false;
		for(int i=1;i<=sum/2;i++)
		{
			for(int j=1;j<=n;j++)
			{
				dp[i][j]=dp[i][j-1];
				if(i>=a[j-1])
					dp[i][j]=dp[i][j-1]||dp[i-a[j-1]][j-1];
			}
		}
		return dp[sum/2][n];
	}	
	int partitionDriver()
	{
		int n;
		cin>>n;
		int a[n],sum=0;
		for(int i=0;i<n;i++)
		{
			cin>>a[i];
			sum+=a[i];
		}	
		cout<<partition(a,n,sum);
	}
	//
