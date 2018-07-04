#ifndef itkPointSetPropertiesCalculator_h
#define itkPointSetPropertiesCalculator_h

#include <itkPointSet.h>
#include <itkNumericTraits.h>

namespace itk
{
template< typename TPointSet >
class PointSetPropertiesCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef PointSetPropertiesCalculator< TPointSet >    Self;
  typedef Object                                       Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetPropertiesCalculator, Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(Dimension, unsigned int, TPointSet::PointDimension);

  /** Standard scalar type within this class. */
  typedef double ScalarType;

  /** Standard types and pointers within this class. */
  typedef TPointSet PointSetType;
  typedef typename PointSetType::ConstPointer                   PointSetConstPointer;
  typedef typename PointSetType::PointType                      PointType;
  typedef typename PointSetType::PointsContainer::ConstPointer  PointsContainerConstPointer;
  typedef typename PointSetType::PointsContainerConstIterator   IteratorType;

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
    m_NumberOfPoints = m_PointSet->GetNumberOfPoints();

    PointsContainerConstPointer points = m_PointSet->GetPoints();
    m_Center.Fill(itk::NumericTraits<typename PointType::ValueType>::ZeroValue());

    // compute center
    for (IteratorType it = points->Begin(); it != points->End(); ++it) 
    {
      PointType point = it.Value();

      for (size_t n = 0; n < Dimension; ++n) 
      {
        m_Center[n] += point[n];
      }
    }

    for (size_t n = 0; n < Dimension; ++n) 
    {
      m_Center[n] /= m_NumberOfPoints;
    }

    // compute radius
    m_Scale = itk::NumericTraits< ScalarType >::ZeroValue();

    for (IteratorType it = points->Begin(); it != points->End(); ++it) 
    {
      PointType point = it.Value();

      for (size_t n = 0; n < Dimension; ++n) 
      {
        m_Scale += pow(point[n] - m_Center[n], 2);
      }
    }

    m_Scale = sqrt(m_Scale / m_NumberOfPoints);
    m_Valid = true;
  }

protected:
  PointSetPropertiesCalculator() {}
  virtual ~PointSetPropertiesCalculator() {};

  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE
  {
    Superclass::PrintSelf(os, indent);
    os << std::endl;
    os << indent << "Point set " << m_PointSet.GetPointer() << std::endl;
    os << indent << "Number    " << m_PointSet->GetNumberOfPoints() << std::endl;
    os << indent << "Center    " << m_Center << std::endl;
    os << indent << "Scale     " << m_Scale << std::endl;
    os << std::endl;
  }

  size_t m_NumberOfPoints;
  PointSetConstPointer m_PointSet;
  PointType m_Center;
  ScalarType m_Scale;
  bool m_Valid = false;

private:
  PointSetPropertiesCalculator(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
}

#endif
