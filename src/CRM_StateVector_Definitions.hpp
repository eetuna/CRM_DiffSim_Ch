#pragma once
#include <assert.h>
//#include <typeinfo>


constexpr unsigned int CURRENT_ACT_VECTOR_DIM = NUM_ACT_SET * 3;
constexpr unsigned int ACT_VECTOR_DIM = CURRENT_ACT_VECTOR_DIM + 1;
constexpr unsigned int TIP_FORCE_DIM = 3;
constexpr unsigned int SIMPLE_STATE_VECTOR_SIZE = (3 + 9 + 3);		// p, R, u
#ifdef ANALYTICAL_SE3_STEP
constexpr unsigned int SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE = SIMPLE_STATE_VECTOR_SIZE - (3 + 9);  // exclude p and R
#else
constexpr unsigned int SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE = SIMPLE_STATE_VECTOR_SIZE;  // include p and R
#endif
//constexpr unsigned int IVP_JACOBIAN_STATE_SIZE = (3 + 3 + 3) * (3 + 3 + 3 + NUM_ACT_SET * 3 + 1 + TIP_FORCE_DIM);  // (_p, ws, _u) x (p0, w0, u0, zc, zl, ftip) 
constexpr unsigned int IVP_JACOBIAN_FULL_STATE_SIZE = (3 + 3 + 3) * (3 + NUM_ACT_SET * 3 + 1 + TIP_FORCE_DIM);  // (_p, ws, _u) x (u0, zc, zl, ftip) 
constexpr unsigned int IVP_JACOBIAN_MINI_STATE_SIZE = (3 + 3 + 3) * (3);  // (_p, ws, _u) x (u0) 


// type definition for mapping raw arrays to Eigen matrices
//   -- note that, by convention, we have been using row major order for matrices
template<typename adType, int Rows, int Cols>
using __EMT = Eigen::Map<Eigen::Matrix<adType, Rows, Cols, Eigen::RowMajor> >;

template<typename adType, int Rows>
using __EVT = Eigen::Map<Eigen::Matrix<adType, Rows, 1> >;



template <typename adType, unsigned int N>
class StateBase {
public:
	StateBase() : _dim(N) {};
	StateBase(const StateBase<adType, N>& obj) : _dim(N) {  // copy constructor - it is a must have if we have dynamic memory allocation
		for (unsigned int i = 0; i < N; i++) data[i] = obj.data[i];
	};
	virtual ~StateBase() {};

	//adType& operator[](const unsigned int idx) { assert(idx < N);  return data[idx]; };
	//const adType& operator[](const unsigned int idx) const { assert(idx < N);  return data[idx]; };
	virtual const StateBase<adType, N> operator+(const StateBase<adType, N>& obj) const {
		StateBase<adType, N> result;
		for (unsigned int i = 0; i < N; i++) result.data[i] = this->data[i] + obj.data[i];
		return result;
	};	
	const StateBase<adType, N> operator*(const adType& scalar) {
		StateBase<adType, N> result;
		for (unsigned int i = 0; i < N; i++) result.data[i] = this->data[i] * scalar;
		return result;
	};
	friend const StateBase<adType, N> operator*(const adType& scalar, const StateBase<adType, N>& obj) {
		StateBase<adType, N> result;
		for (unsigned int i = 0; i < N; i++) result.data[i] = scalar * obj.data[i];
		return result;
	};

	adType absmax() {
		adType mx = 0.0;
		for (unsigned int i = 0; i < _dim; i++) mx = MAX(mx, abs(data[i]));
		return mx;
	};

	friend std::ostream& operator<<(std::ostream& os, const StateBase<adType, N>& obj) {
		for (unsigned int i = 0; i < obj._dim-1; i++) os << obj.data[i] << ", ";
		os << obj.data[obj._dim - 1] << "\n";
		return os;
	};

protected:
	adType data[N] = { };		// component array, initialized to zero
	unsigned int _dim;

private:

};

#ifdef ANALYTICAL_SE3_STEP
#define __INIT() _u(this->data)
#else
#define __INIT() _p(this->data), _R(this->data+3), _u(this->data+3+9)
#endif

template <typename adType>
class StateDerivativeVector : public StateBase<adType, SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE> {
public:
	StateDerivativeVector() : StateBase<adType, SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE>(), __INIT() { };
	StateDerivativeVector(const StateDerivativeVector<adType>& obj) :	// copy constructor
		StateBase<adType, SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE>(obj), __INIT() { };
	StateDerivativeVector(const StateBase<adType, SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE>& obj) :
		StateBase<adType, SIMPLE_STATE_DERIVATIVE_VECTOR_SIZE>(obj), __INIT() { };
	virtual ~StateDerivativeVector() {};

	virtual StateDerivativeVector<adType>& operator=(const StateDerivativeVector<adType>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		for (unsigned int i = 0; i < this->_dim; i++) this->data[i] = rhs.data[i];
		return *this;
	};

#ifndef ANALYTICAL_SE3_STEP
	adType* const _p;
	adType* const _R;
#endif
	adType* const _u;

protected:

private:

};

#undef __INIT


#define __INIT() _p(this->data), _R(this->data + 3), _u(this->data + 3 + 9)

template <typename adType>
class StateVector : public StateBase<adType, SIMPLE_STATE_VECTOR_SIZE> {
public:
	StateVector() : StateBase<adType, SIMPLE_STATE_VECTOR_SIZE>(), __INIT() { };
	StateVector(const StateVector<adType>& obj) :						// copy constructor
		StateBase<adType, SIMPLE_STATE_VECTOR_SIZE>(obj), __INIT() { };
	StateVector(const StateBase<adType, SIMPLE_STATE_VECTOR_SIZE>& obj) :
		StateBase<adType, SIMPLE_STATE_VECTOR_SIZE>(obj), __INIT() { };
	virtual ~StateVector() {};

	virtual StateVector<adType>& operator=(const StateVector<adType>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		for (unsigned int i = 0; i < this->_dim; i++) this->data[i] = rhs.data[i];
		return *this;
	};
	virtual StateVector<adType>& operator=(const StateDerivativeVector<adType>& rhs) {
		for (unsigned int i = 0; i < 3; i++) _u[i] = rhs._u[i];
		return *this;
	};
	virtual const StateVector<adType> operator+(const StateVector<adType>& obj) const {
		StateVector<adType> result;
		for (unsigned int i = 0; i < this->_dim; i++) result.data[i] = this->data[i] + obj.data[i];
		return result;
	};
	virtual const StateVector<adType> operator+(const StateDerivativeVector<adType>& obj) const {
		StateVector<adType> result;
		for (unsigned int i = 0; i < this->_dim; i++) result.data[i] = this->data[i];
		for (unsigned int i = 0; i < 3; i++) result._u[i] += obj._u[i];
		return result;
	};


	adType* const _p;
	adType* const _R;
	adType* const _u;

protected:

private:

};

#undef __INIT


#define __INIT() _p_u0(this->data+0), _ws_u0(this->data+9), _u_u0(this->data+18)

template <typename adType>
class IVPJacobiansMini : public StateBase<adType, IVP_JACOBIAN_MINI_STATE_SIZE> {
public:
	IVPJacobiansMini() : StateBase<adType, IVP_JACOBIAN_MINI_STATE_SIZE>(), __INIT() { };
	IVPJacobiansMini(const IVPJacobiansMini<adType>& obj) :						// copy constructor
		StateBase<adType, IVP_JACOBIAN_MINI_STATE_SIZE>(obj), __INIT() { };
	IVPJacobiansMini(const StateBase<adType, IVP_JACOBIAN_MINI_STATE_SIZE>& obj) :
		StateBase<adType, IVP_JACOBIAN_MINI_STATE_SIZE>(obj), __INIT() { };
	virtual ~IVPJacobiansMini() {};

	virtual IVPJacobiansMini<adType>& operator=(const IVPJacobiansMini<adType>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		for (unsigned int i = 0; i < this->_dim; i++) this->data[i] = rhs.data[i];
		return *this;
	};

	friend std::ostream& operator<<(std::ostream& os, const IVPJacobiansMini<adType>& obj) {
		os << "JIVP_p_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_p_u0(obj._p_u0);
		os << JIVP_p_u0 << "];\n";
		os << "JIVP_ws_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_ws_u0(obj._ws_u0);
		os << JIVP_ws_u0 << "];\n";
		os << "JIVP_u_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_u_u0(obj._u_u0);
		os << JIVP_u_u0 << "];\n";

		return os;
	};

	adType* const _p_u0;
	adType* const _ws_u0;
	adType* const _u_u0;

protected:

private:

};

#undef __INIT



//#define __INIT() _p_p0(this->data + 0), _p_w0(this->data + 9), _p_u0(this->data + 18), _ws_p0(this->data + 27), _ws_w0(this->data + 36), _ws_u0(this->data + 45), _u_p0(this->data + 54), _u_w0(this->data + 63), _u_u0(this->data + 72), _p_zl(this->data + 81), _ws_zl(this->data + 84), _u_zl(this->data + 87), _p_ft(this->data + 90), _ws_ft(this->data + 99), _u_ft(this->data + 108), _p_zc(this->data + 117), _ws_zc(this->data + 117 + 3 * CURRENT_ACT_VECTOR_DIM), _u_zc(this->data + 117 + 6 * CURRENT_ACT_VECTOR_DIM)
#define __INIT() _p_u0(this->data+0), _ws_u0(this->data+9), _u_u0(this->data+18), _p_zl(this->data+27), _ws_zl(this->data+30), _u_zl(this->data+33), _p_ft(this->data+36), _ws_ft(this->data+45), _u_ft(this->data+54), _p_zc(this->data+63), _ws_zc(this->data+63+3*CURRENT_ACT_VECTOR_DIM), _u_zc(this->data+63+6*CURRENT_ACT_VECTOR_DIM)

template <typename adType>
class IVPJacobiansFull : public StateBase<adType, IVP_JACOBIAN_FULL_STATE_SIZE> {
public:
	IVPJacobiansFull() : StateBase<adType, IVP_JACOBIAN_FULL_STATE_SIZE>(), __INIT() { };
	IVPJacobiansFull(const IVPJacobiansFull<adType>& obj) :						// copy constructor
		StateBase<adType, IVP_JACOBIAN_FULL_STATE_SIZE>(obj), __INIT() { };
	IVPJacobiansFull(const StateBase<adType, IVP_JACOBIAN_FULL_STATE_SIZE>& obj) :
		StateBase<adType, IVP_JACOBIAN_FULL_STATE_SIZE>(obj), __INIT() { };
	virtual ~IVPJacobiansFull() {};

	virtual IVPJacobiansFull<adType>& operator=(const IVPJacobiansFull<adType>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		for (unsigned int i = 0; i < this->_dim; i++) this->data[i] = rhs.data[i];
		return *this;
	};

	friend std::ostream& operator<<(std::ostream& os, const IVPJacobiansFull<adType>& obj) {
		os << "JIVP_p_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_p_u0(obj._p_u0);
		os << JIVP_p_u0 << "];\n";
		os << "JIVP_p_zc=[ \n";
		__EMT<adType, 3, 3 * NUM_ACT_SET> JIVP_p_zc(obj._p_zc);
		os << JIVP_p_zc << "];\n";
		os << "JIVP_p_zl=[ \n";
		__EVT<adType, 3> JIVP_p_zl(obj._p_zl);
		os << JIVP_p_zl << "];\n";
		os << "JIVP_p_ft=[ \n";
		__EMT<adType, 3, 3> JIVP_p_ft(obj._p_ft);
		os << JIVP_p_ft << "];\n";
		os << "\n";
		os << "JIVP_ws_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_ws_u0(obj._ws_u0);
		os << JIVP_ws_u0 << "];\n";
		os << "JIVP_ws_zc=[ \n";
		__EMT<adType, 3, 3 * NUM_ACT_SET> JIVP_ws_zc(obj._ws_zc);
		os << JIVP_ws_zc << "];\n";
		os << "JIVP_ws_zl=[ \n";
		__EVT<adType, 3> JIVP_ws_zl(obj._ws_zl);
		os << JIVP_ws_zl << "];\n";
		os << "JIVP_ws_ft=[ \n";
		__EMT<adType, 3, 3> JIVP_ws_ft(obj._ws_ft);
		os << JIVP_ws_ft << "];\n";
		os << "\n";
		os << "JIVP_u_u0=[ \n";
		__EMT<adType, 3, 3> JIVP_u_u0(obj._u_u0);
		os << JIVP_u_u0 << "];\n";
		os << "JIVP_u_zc=[ \n";
		__EMT<adType, 3, 3 * NUM_ACT_SET> JIVP_u_zc(obj._u_zc);
		os << JIVP_u_zc << "];\n";
		os << "JIVP_u_zl=[ \n";
		__EVT<adType, 3> JIVP_u_zl(obj._u_zl);
		os << JIVP_u_zl << "];\n";
		os << "JIVP_u_ft=[ \n";
		__EMT<adType, 3, 3> JIVP_u_ft(obj._u_ft);
		os << JIVP_u_ft << "];\n";

		return os;
	};

	adType* const _p_u0;
	adType* const _p_zc;
	adType* const _p_zl;
	adType* const _p_ft;
	adType* const _ws_u0;
	adType* const _ws_zc;
	adType* const _ws_zl;
	adType* const _ws_ft;
	adType* const _u_u0;
	adType* const _u_zc;
	adType* const _u_zl;
	adType* const _u_ft;

protected:

private:

};

#undef __INIT



template <typename adType, template<typename> typename IVPJacobians>
class AugmentedStateDerivativeVector : public StateDerivativeVector<adType>, public IVPJacobians<adType> {
public:
	AugmentedStateDerivativeVector() : StateDerivativeVector<adType>(), IVPJacobians<adType>() { };
	AugmentedStateDerivativeVector(const AugmentedStateDerivativeVector<adType, IVPJacobians>& obj) :						// copy constructor
		StateDerivativeVector<adType>(obj), IVPJacobians<adType>(obj) { };
	virtual ~AugmentedStateDerivativeVector() {};

	virtual AugmentedStateDerivativeVector<adType, IVPJacobians>& operator=(const AugmentedStateDerivativeVector<adType, IVPJacobians>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		StateDerivativeVector<adType>::operator=(rhs);
		IVPJacobians<adType>::operator=(rhs);
		return *this;
	};
	virtual const AugmentedStateDerivativeVector<adType, IVPJacobians> operator+(const AugmentedStateDerivativeVector<adType, IVPJacobians>& obj) const {
		AugmentedStateDerivativeVector<adType, IVPJacobians> result;
		static_cast<StateDerivativeVector<adType>&>(result) = StateDerivativeVector<adType>::operator+(obj);
		static_cast<IVPJacobians<adType>&>(result) = IVPJacobians<adType>::operator+(obj);
		return result;
	};
	const AugmentedStateDerivativeVector<adType, IVPJacobians> operator*(const adType& scalar) {
		AugmentedStateDerivativeVector<adType, IVPJacobians> result;
		static_cast<StateDerivativeVector<adType>&>(result) = StateDerivativeVector<adType>::operator*(scalar);
		static_cast<IVPJacobians<adType>&>(result) = IVPJacobians<adType>::operator*(scalar);
		return result;
	};
	friend const AugmentedStateDerivativeVector<adType, IVPJacobians> operator*(const adType& scalar, const AugmentedStateDerivativeVector<adType, IVPJacobians>& obj) {
		AugmentedStateDerivativeVector<adType, IVPJacobians> result;
		static_cast<StateDerivativeVector<adType>&>(result) = scalar * static_cast<const StateDerivativeVector<adType>&>(obj);
		static_cast<IVPJacobians<adType>&>(result) = scalar * static_cast<const IVPJacobians<adType>&>(obj);
		return result;
	};

	friend std::ostream& operator<<(std::ostream& os, const AugmentedStateDerivativeVector<adType, IVPJacobians>& obj) {
		os << static_cast<const StateDerivativeVector<adType>&>(obj) << "\n" << static_cast<const IVPJacobians<adType>&>(obj);
		return os;
	};


protected:

private:

};


template <typename adType, template<typename> typename IVPJacobians>
class AugmentedStateVector : public StateVector<adType>, public IVPJacobians<adType> {
public:
	AugmentedStateVector() : StateVector<adType>(), IVPJacobians<adType>() { };
	AugmentedStateVector(const AugmentedStateVector<adType, IVPJacobians>& obj) :						// copy constructor
		StateVector<adType>(obj), IVPJacobians<adType>(obj) { };
	virtual ~AugmentedStateVector() {};


	virtual AugmentedStateVector<adType, IVPJacobians>& operator=(const AugmentedStateVector<adType, IVPJacobians>& rhs) {
		if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.
		StateVector<adType>::operator=(rhs);
		IVPJacobians<adType>::operator=(rhs);
		return *this;
	};
	virtual const AugmentedStateVector<adType, IVPJacobians> operator+(const AugmentedStateVector<adType, IVPJacobians>& obj) const {
		AugmentedStateVector<adType, IVPJacobians> result;
		static_cast<StateVector<adType>&>(result) = StateVector<adType>::operator+(static_cast<const StateVector<adType>&>(obj));
		static_cast<IVPJacobians<adType>&>(result) = IVPJacobians<adType>::operator+(static_cast<const IVPJacobians<adType>&>(obj));
		return result;
	};
	virtual const AugmentedStateVector<adType, IVPJacobians> operator*(const adType& scalar) {
		AugmentedStateVector<adType, IVPJacobians> result;
		static_cast<StateVector<adType>&>(result) = StateVector<adType>::operator*(scalar);
		static_cast<IVPJacobians<adType>&>(result) = IVPJacobians<adType>::operator*(scalar);
		return result;
	};
	friend const AugmentedStateVector<adType, IVPJacobians> operator*(const adType& scalar, const AugmentedStateVector<adType, IVPJacobians>& obj) {
		AugmentedStateVector<adType, IVPJacobians> result;
		static_cast<StateVector<adType>&>(result) = scalar * static_cast<const StateVector<adType>&>(obj);
		static_cast<IVPJacobians<adType>&>(result) = scalar * static_cast<const IVPJacobians<adType>&>(obj);
		return result;
	};
	virtual AugmentedStateVector<adType, IVPJacobians>& operator=(const AugmentedStateDerivativeVector<adType, IVPJacobians>& rhs) {
		//if (this == &rhs)  return *this;		// If same object skip assignment, and just return *this.  // they cannot be the same, so no need for this
		StateVector<adType>::operator=(static_cast<const StateDerivativeVector<adType>&>(rhs));
		IVPJacobians<adType>::operator=(rhs);
		return *this;
	};
	virtual const AugmentedStateVector<adType, IVPJacobians> operator+(const AugmentedStateDerivativeVector<adType, IVPJacobians>& obj) const {
		AugmentedStateVector<adType, IVPJacobians> result;
		static_cast<StateVector<adType>&>(result) = StateVector<adType>::operator+(static_cast<const StateDerivativeVector<adType>&>(obj));
		static_cast<IVPJacobians<adType>&>(result) = IVPJacobians<adType>::operator+(obj);
		return result;
	};

	adType absmax() {
		return MAX(static_cast<StateVector<adType>*>(this)->absmax(), static_cast<IVPJacobians<adType>*>(this)->absmax());
	};


	friend std::ostream& operator<<(std::ostream& os, const AugmentedStateVector<adType, IVPJacobians>& obj) {
		os << static_cast<const StateVector<adType>&>(obj) << "\n" << static_cast<const IVPJacobians<adType>&>(obj);
		return os;
	};


protected:

private:

};


