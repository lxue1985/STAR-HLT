// data.h: interface for the data class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATA_H__FAB4AC21_C059_402C_B493_845F7B4D3382__INCLUDED_)
#define AFX_DATA_H__FAB4AC21_C059_402C_B493_845F7B4D3382__INCLUDED_
class data  

{

    public:
    int midlength;
    int put(char * msource,int mlength);
    int put(char * msource1,int mlength1,char * msource2,int mlength2,char * msource3,int mlength3,char * msource4,int mlength4,int mtype);
    int put(char * msource1,int mlength1,char * msource2,int mlength2,char * msource3,int mlength3,int mtype);
    int put(char *,int,char *,int ,int);
    int get();
    int put(char * msource,int mlength,int mtype);
    void end();
    void start();
    int length;
    int mlength;
    int number;
    char * result;
    char * source;
    int source_length;
    int type;
    data();
    data(int aa);
    virtual ~data();
    int max_length;
    int max_number;
    
}
;

#endif // !defined(AFX_DATA_H__FAB4AC21_C059_402C_B493_845F7B4D3382__INCLUDED_)

