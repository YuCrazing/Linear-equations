
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

    int GuassInOrder(); //˳���˹��Ԫ��
    int GuassColumnPrincipleComponent();    //����Ԫ��˹��Ԫ��
    int SquareRoot();   //ƽ������
    int Chasing();  //׷�Ϸ�
};
