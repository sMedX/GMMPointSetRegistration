/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkGMMPointSetToPointSetMetricBase_h
#define itkGMMPointSetToPointSetMetricBase_h

#include "itkSingleValuedCostFunction.h"
#include "itkTransform.h"
#include "itkPointSet.h"
#include "itkMacro.h"
#include "itkPointsLocator.h"

namespace itk
{
/** \class PointSetToPointSetMetric
 * \brief Computes similarity between two point sets.
 *
 * This Class is templated over the type of the two point-sets.  It
 * expects a Transform to be plugged in.  This particular
 * class is the base class for a hierarchy of point-set to point-set metrics.
 *
 * This class computes a value that measures the similarity between the fixed point-set
 * and the transformed moving point-set.
 *
 * \ingroup RegistrationMetrics
 *
 * \ingroup ITKRegistrationCommon
 */

template< typename TFixedPointSet,  typename TMovingPointSet >
class GMMPointSetToPointSetMetricBase :public SingleValuedCostFunction
{
public:

  /** Standard class typedefs. */
  typedef GMMPointSetToPointSetMetricBase   Self;
  typedef SingleValuedCostFunction          Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Type used for representing point components  */
  typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMPointSetToPointSetMetricBase, SingleValuedCostFunction);

  /**  Type of the moving Pointset. */
  typedef TMovingPointSet                           MovingPointSetType;
  typedef typename MovingPointSetType::PointType    MovingPointType;
  typedef typename MovingPointSetType::ConstPointer MovingPointSetConstPointer;

  /**  Type of the fixed Pointset. */
  typedef TFixedPointSet                            FixedPointSetType;
  typedef typename FixedPointSetType::PointType     FixedPointType;
  typedef typename FixedPointSetType::ConstPointer  FixedPointSetConstPointer;

  /** Constants for the pointset dimensions */
  itkStaticConstMacro(MovingPointSetDimension, unsigned int, TMovingPointSet::PointDimension);
  itkStaticConstMacro(FixedPointSetDimension, unsigned int, TFixedPointSet::PointDimension);
  itkStaticConstMacro(PointDimension, unsigned int, TMovingPointSet::PointDimension);

  typedef typename FixedPointSetType::PointsContainer                     FixedPointsContainer;
  typedef typename FixedPointSetType::PointsContainer::Pointer            FixedPointsPointer;
  typedef typename FixedPointsContainer::ConstIterator                    FixedPointIterator;
  typedef typename FixedPointSetType::PointDataContainer::ConstIterator   FixedPointDataIterator;
  typedef itk::PointsLocator<FixedPointsContainer>                        FixedPointsLocatorType;
  typedef typename FixedPointsLocatorType::NeighborsIdentifierType        FixedNeighborsIdentifierType;
  typedef typename FixedNeighborsIdentifierType::const_iterator           FixedNeighborsIteratorType;

  typedef typename MovingPointSetType::PointsContainer                    MovingPointsContainer;
  typedef typename MovingPointSetType::PointsContainer::Pointer           MovingPointsPointer;
  typedef typename MovingPointsContainer::ConstIterator                   MovingPointIterator;
  typedef typename MovingPointSetType::PointDataContainer::ConstIterator  MovingPointDataIterator;
  typedef itk::PointsLocator<MovingPointsContainer>                       MovingPointsLocatorType;
  typedef typename MovingPointsLocatorType::NeighborsIdentifierType       MovingNeighborsIdentifierType;
  typedef typename MovingNeighborsIdentifierType::const_iterator          MovingNeighborsIteratorType;

  /**  Type of the Transform Base class */
  typedef Transform< CoordinateRepresentationType,
                     itkGetStaticConstMacro(MovingPointSetDimension),
                     itkGetStaticConstMacro(FixedPointSetDimension) > TransformType;

  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType                       DerivativeType;
  typedef DerivativeType::ValueType                        DerivativeValueType;
  typedef FixedArray<DerivativeValueType, PointDimension>  LocalDerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  /** Get/Set the scale.  */
  itkSetMacro(Scale, double);
  itkGetMacro(Scale, double);

  /** Get/Set the Fixed point set.  */
  itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
  itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

  /** Get/Set the Moving point set.  */
  itkSetConstObjectMacro(MovingPointSet, MovingPointSetType);
  itkGetConstObjectMacro(MovingPointSet, MovingPointSetType);

  /** Get/Set boolean flag to initialize KdTree.  */
  itkSetMacro(UseFixedPointSetKdTree, bool);
  itkGetMacro(UseFixedPointSetKdTree, bool);

  itkSetMacro(UseMovingPointSetKdTree, bool);
  itkGetMacro(UseMovingPointSetKdTree, bool);

  itkSetMacro(Radius, double);
  itkGetMacro(Radius, double);

  /** Connect the Transform. */
  itkSetObjectMacro(Transform, TransformType);

  /** Get a pointer to the Transform.  */
  itkGetModifiableObjectMacro(Transform, TransformType);

  /** Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters, DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & Value, DerivativeType & Derivative) const ITK_OVERRIDE;

  /** Calculates the local metric value for a single point. */
  virtual MeasureType GetLocalNeighborhoodValue(const MovingPointType & point) const = 0;

  /** Calculates the local value/derivative for a single point.*/
  virtual void GetLocalNeighborhoodValueAndDerivative(const MovingPointType &, MeasureType &, LocalDerivativeType &) const = 0;

  /** Initialize to prepare for a particular iteration, generally an iteration of optimization. Distinct from Initialize()
  * which is a one-time initialization. */
  virtual void InitializeForIteration(const ParametersType & parameters) const;

  /** Set the parameters defining the Transform. */
  void SetTransformParameters(const ParametersType & parameters) const;

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
  { return m_Transform->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void)
  throw ( ExceptionObject );

protected:
  GMMPointSetToPointSetMetricBase();
  virtual ~GMMPointSetToPointSetMetricBase() {}
  void InitializeFixedTree();
  void InitializeMovingTree();
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  FixedPointSetConstPointer m_FixedPointSet;
  MovingPointSetConstPointer m_MovingPointSet;
  mutable typename MovingPointSetType::Pointer m_TransformedMovingPointSet;

  mutable TransformPointer m_Transform;
  size_t m_NumberOfParameters;

  mutable TransformJacobianType m_Jacobian;
  mutable TransformJacobianType m_JacobianCache;

  double m_Scale;

  size_t m_NumberOfFixedPoints;
  size_t m_NumberOfMovingPoints;

  double m_NormalizingValueFactor;
  double m_NormalizingDerivativeFactor;

  typename FixedPointsLocatorType::Pointer   m_FixedPointsLocator;
  typename MovingPointsLocatorType::Pointer  m_MovingPointsLocator;
  bool m_UseFixedPointSetKdTree;
  bool m_UseMovingPointSetKdTree;
  double m_Radius;

private:
  GMMPointSetToPointSetMetricBase(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMPointSetToPointSetMetricBase.hxx"
#endif

#endif
