#ifndef itkGMMMLEPointSetToPointSetMetric_h
#define itkGMMMLEPointSetToPointSetMetric_h

#include <itkCovariantVector.h>
#include <itkPoint.h>
#include "itkGMMPointSetToPointSetMetricBase.h"

namespace itk
{
/** \class GMMRigidPointSetToPointSetMetric
 * \brief Computes similarity between pixel values of a point set and
 * intensity values of an image.
 *
 * This metric computes the average squared differences between pixels
 * in the point set and the transformed set of pixels.
 *
 * Spatial correspondence between both images is established through a
 * Transform.
 */
template< typename TFixedPointSet, typename TMovingPointSet >
class GMMMLEPointSetToPointSetMetric : public GMMPointSetToPointSetMetricBase < TFixedPointSet, TMovingPointSet >
{
public:
  /** Standard class typedefs. */
  typedef GMMMLEPointSetToPointSetMetric                                     Self;
  typedef GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet > Superclass;
  typedef SmartPointer< Self >                                               Pointer;
  typedef SmartPointer< const Self >                                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMMLEPointSetToPointSetMetric, GMMPointSetToPointSetMetricBase);

  /** Types transferred from the base class */
  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::TransformParametersType    TransformParametersType;
  typedef typename Superclass::TransformJacobianType      TransformJacobianType;
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::FixedPointSetType          FixedPointSetType;
  typedef typename Superclass::MovingPointSetType         MovingPointSetType;
  typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
  typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;
  typedef typename Superclass::FixedPointIterator         FixedPointIterator;
  typedef typename Superclass::MovingPointIterator        MovingPointIterator;
  typedef typename Superclass::LocalDerivativeType        LocalDerivativeType;

  /** Initialize the Metric by making sure that all the components are present and plugged together correctly.*/
  virtual void Initialize() throw (ExceptionObject)ITK_OVERRIDE;

protected:
  GMMMLEPointSetToPointSetMetric();
  virtual ~GMMMLEPointSetToPointSetMetric() {}

  /** Calculates the local metric value for a single point.*/
  virtual MeasureType GetLocalNeighborhoodValue(const MovingPointIterator &) const ITK_OVERRIDE;

  /** Calculates the local value/derivative for a single point.*/
  virtual bool GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator &, MeasureType &, LocalDerivativeType &) const ITK_OVERRIDE;

private:
  GMMMLEPointSetToPointSetMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMMLEPointSetToPointSetMetric.hxx"
#endif

#endif
