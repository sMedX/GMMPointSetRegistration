#ifndef itkNormalizePointSet_h
#define itkNormalizePointSet_h

#include <itkPointSet.h>
#include <itkNumericTraits.h>
#include <itkScaleTransform.h>
#include "itkPointSetPropertiesCalculator.h"

namespace itk
{
template< typename TPointSet >
class NormalizePointSet : public Object
{
public:
  /** Standard class typedefs. */
  typedef NormalizePointSet< TPointSet >    Self;
  typedef Object                            Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizePointSet, Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(Dimension, unsigned int, TPointSet::PointDimension);

  /** Standard scalar type within this class. */
  typedef double ScalarType;

  /** Standard types and pointers within this class. */
  typedef TPointSet PointSetType;
  typedef typename PointSetType::Pointer                        PointSetPointer;
  typedef typename PointSetType::ConstPointer                   PointSetConstPointer;
  typedef typename PointSetType::PointType                      PointType;
  typedef typename PointSetType::PointsContainer                PointsContainer;
  typedef typename PointSetType::PointsContainer::Pointer       PointsContainerPointer;
  typedef typename PointSetType::PointsContainer::ConstPointer  PointsContainerConstPointer;
  typedef typename PointSetType::PointsContainerConstIterator   IteratorType;
  typedef itk::ScaleTransform<double, Dimension> TransformType;

  typedef PointSetPropertiesCalculator<PointSetType>            PointSetPropertiesCalculatorType;

  itkGetObjectMacro(Transform, TransformType);

  /** Set the input image. */
  virtual void SetPointSet(const PointSetType *points)
  {
    if ( m_PointSet != points ) 
    {
      m_PointSet = points;
      this->Modified();
      m_Valid = false;
    }
  }

  /** Get output point set.*/
  PointSetPointer GetOutput() const
  {
    if (!m_Valid) {
      itkExceptionMacro(<< "GetOutput() invoked, but the properties have not been computed. Call Compute() first.");
    }
    return m_OutputPointSet;
  }

  /** Get scale.*/
  ScalarType GetScale() const
  {
    if (!m_Valid) {
      itkExceptionMacro(<< "GetScale() invoked, but the properties have not been computed. Call Compute() first.");
    }
    return m_Scale;
  }

  /** Get center in physical coordinates.*/
  PointType GetCenter() const
  {
    if (!m_Valid) {
      itkExceptionMacro(<< "GetCenter() invoked, but the properties have not been computed. Call Compute() first.");
    }
    return m_Center;
  }

  void Compute()
  {
    m_PointSetCalculator = PointSetPropertiesCalculatorType::New();
    m_PointSetCalculator->SetPointSet(m_PointSet);
    m_PointSetCalculator->Compute();
    m_Center = m_PointSetCalculator->GetCenter();
    m_Scale = m_PointSetCalculator->GetScale();

    m_Transform = TransformType::New();
    m_Transform->SetScale(1.0 / m_Scale);
    m_Transform->SetCenter(m_Center);

    typename TransformType::TranslationType translation;
    for (unsigned int n = 0; n < Dimension; ++n) {
      translation[n] = -m_Center[n];
    }

    m_Transform->SetTranslation(translation);

    typename PointsContainer::Pointer points = PointsContainer::New();
    m_OutputPointSet = PointSetType::New();

    // normalize point set
    for (IteratorType it = m_PointSet->GetPoints()->Begin(); it != m_PointSet->GetPoints()->End(); ++it)
    {
      points->InsertElement(it.Index(), m_Transform->TransformPoint(it.Value()));
    }

    m_OutputPointSet->SetPoints(points);
    m_Valid = true;
  }

  void PrintReport(std::ostream& os)
  {
    os << "points " << m_PointSet->GetNumberOfPoints() << std::endl;
    os << "center " << m_Center << std::endl;
    os << "scale  " << m_Scale << std::endl;
    os << std::endl;
  }

protected:
  NormalizePointSet() {}
  virtual ~NormalizePointSet() {};

  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "PointSet: " << m_PointSet.GetPointer() << std::endl;
    os << indent << "PointSet: " << m_OutputPointSet.GetPointer() << std::endl;
  }

  typename PointSetPropertiesCalculatorType::Pointer m_PointSetCalculator;
  PointSetConstPointer m_PointSet;
  PointSetPointer m_OutputPointSet;
  PointsContainerPointer points;
  bool m_Valid = false;

  PointType m_Center;
  ScalarType m_Scale;
  typename TransformType::Pointer m_Transform;

private:
  NormalizePointSet(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
}

#endif
