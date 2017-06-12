#ifndef __agtkObservers_h
#define __agtkObservers_h

#include <itkCommand.h>
#include <itkLBFGSBOptimizer.h>

namespace agtk
{
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() {};
public:
  typedef itk::LBFGSBOptimizer OptimizerType;
  typedef const OptimizerType *                            OptimizerPointer;
  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
    Execute((const itk::Object *)caller, event);
  }
  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
    OptimizerPointer optimizer = static_cast<OptimizerPointer>(object);
    if (!itk::IterationEvent().CheckEvent(&event)) {
      return;
    }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};
}
#endif // __agtkObservers_h
