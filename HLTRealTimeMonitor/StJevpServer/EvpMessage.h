#ifndef _EVPMESSAGE_H_
#define _EVPMESSAGE_H_

#include <TROOT.h>

class EvpMessage : public TObject {
 public:
  char *cmd;
  char *args;

  void setCmd(char *cmd) {
    if(this->cmd) delete this->cmd;
    int len = strlen(cmd);
    this->cmd = new char[len+1];
    strcpy(this->cmd, cmd);
  };

  void setArgs(char *args) {
    if(this->args) delete this->args;
    int len = strlen(args);
    this->args = new char[len+1];
    strcpy(this->args, args);
  };
  

  EvpMessage() {
    cmd = NULL;
    args = NULL;
  }

  ~EvpMessage() {
    if(cmd != NULL) delete cmd;
    if(args != NULL) delete args;
  }
  
  ClassDef(EvpMessage, 1);
};

#endif
