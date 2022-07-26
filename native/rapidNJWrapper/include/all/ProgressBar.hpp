#pragma once
#include <stack>
#include <utility>
#include "stdinclude.h"

typedef void (*progress_callback)(double);

class ProgressBar
{
public:
  ProgressBar(void);
  ProgressBar(progress_callback progressCallback);
  void childProgress(double fraction);
  void setProgress(double progress);
  void finish();
  ~ProgressBar(void);
  
private:  
  std::stack<std::pair<double,double> > intervals;
  double curStart;
  double curEnd;
  double currentProgress;
  progress_callback callback;
  bool hasCallback;
};

