
class Matrix{
    int N;
    double a[100][100];
    double b[100];
    double x[100];

public:
    void init();
    void showN();
    void showA();
    void showb();
    void showx();
    void show();

    void swapRow(int _x, int _y);

    int GuassInOrder(); //顺序高斯消元法
    int GuassColumnPrincipleComponent();    //列主元高斯消元法
    int SquareRoot();   //平方根法
    int Chasing();  //追赶法
};
