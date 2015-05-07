#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RUNGE_KUTTA4_IMPLICIT_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_RUNGE_KUTTA4_IMPLICIT_HPP_INCLUDED

#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<
class State ,
class Value = double ,
class Deriv = State ,
class Time = Value ,
class Algebra = typename algebra_dispatcher< State >::algebra_type ,
class Operations = typename operations_dispatcher< State >::operations_type ,
class Resizer = initially_resizer
>

class runge_kutta4_implicit
: public explicit_stepper_base<
  runge_kutta4_implicit< State , Value , Deriv , Time , Algebra , Operations , Resizer > ,
  4 , State , Value , Deriv , Time , Algebra , Operations , Resizer >
{
public:

    typedef explicit_stepper_base<
    runge_kutta4_implicit< State , Value , Deriv , Time , Algebra , Operations , Resizer > ,
    4 , State , Value , Deriv , Time , Algebra , Operations , Resizer > stepper_base_type;

    typedef typename stepper_base_type::state_type state_type;
    typedef typename stepper_base_type::value_type value_type;
    typedef typename stepper_base_type::deriv_type deriv_type;
    typedef typename stepper_base_type::time_type time_type;
    typedef typename stepper_base_type::algebra_type algebra_type;
    typedef typename stepper_base_type::operations_type operations_type;
    typedef typename stepper_base_type::resizer_type resizer_type;

    #ifndef DOXYGEN_SKIP
    typedef typename stepper_base_type::stepper_type stepper_type;
    typedef typename stepper_base_type::wrapped_state_type wrapped_state_type;
    typedef typename stepper_base_type::wrapped_deriv_type wrapped_deriv_type;
    #endif // DOXYGEN_SKIP


    runge_kutta4_implicit( const algebra_type &algebra = algebra_type() ) : stepper_base_type( algebra )
    { }


    template< class System , class StateIn , class DerivIn , class StateOut >
	void do_step_impl( System system , const StateIn &in , DerivIn &dxdt , time_type t , StateOut &out , time_type dt )
	{
        static const value_type val1 = static_cast< value_type >( 1.0 );
		//std::cout << "rk4imp called" << std::endl;

		// use the trivial start values k1 = k2 = f( t_old , x_old )
		state_type k1 = dxdt, k2 = dxdt;
		// perform three fixed point iterations
		iterate_functional( system, k1, k2, dt, t , in );
		iterate_functional( system, k1, k2, dt, t , in );
		iterate_functional( system, k1, k2, dt, t , in );
		// out = in + h/2*(k1 + k2)
		stepper_base_type::m_algebra.for_each4( out, in, k1, k2, typename operations_type::template scale_sum3< value_type , time_type, time_type >( val1, dt/2.0 , dt/2.0 ) );
	}

	template< class System >
	void iterate_functional( System& system, state_type& k1, state_type&k2, time_type dt, time_type t, const state_type& x )
	{
        static const value_type val1 = static_cast< value_type >( 1 );
		static const value_type c1  = static_cast< value_type >( 1.0/2.0 - std::sqrt(3)/6.0 );
		static const value_type c2  = static_cast< value_type >( 1.0/2.0 + std::sqrt(3)/6.0 );
		static const value_type a21 = static_cast< value_type >( 1.0/4.0 - std::sqrt(3)/6.0 );
		static const value_type a12 = static_cast< value_type >( 1.0/4.0 + std::sqrt(3)/6.0 );
		static const value_type a11 = static_cast< value_type >( 1.0/4.0 );
		static const value_type a22 = static_cast< value_type >( 1.0/4.0 );
		// calculate m_x_tmp = x+h(a11*k1+a12*k2)
        stepper_base_type::m_algebra.for_each4( m_x_tmp.m_v , x, k1 , k2 ,
												typename operations_type::template scale_sum3< value_type, time_type , time_type >( val1, a11*dt , a12*dt ) );
		// update k1 = f( t + c1*h, m_x_tmp )
        typename odeint::unwrap_reference< System >::type &sys = system;
		sys( m_x_tmp.m_v , k1 , t+dt*c1 );
		// calculate m_x_tmp = x+h(a21*k1+a22*k2)
        stepper_base_type::m_algebra.for_each4( m_x_tmp.m_v , x, k1 , k2 ,
												typename operations_type::template scale_sum3< value_type, time_type , time_type >( val1, a21*dt , a22*dt ) );
		// update k2 = f( t + c1*h, m_x_tmp )
		sys( m_x_tmp.m_v , k2 , t+dt*c2 );
	}

    template< class StateType >
    void adjust_size( const StateType &x )
    {
        resize_impl( x );
        stepper_base_type::adjust_size( x );
    }

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        bool resized = false;
        resized |= adjust_size_by_resizeability( m_x_tmp , x , typename is_resizeable<state_type>::type() );
        resized |= adjust_size_by_resizeability( m_dxm , x , typename is_resizeable<deriv_type>::type() );
        resized |= adjust_size_by_resizeability( m_dxt , x , typename is_resizeable<deriv_type>::type() );
        resized |= adjust_size_by_resizeability( m_dxh , x , typename is_resizeable<deriv_type>::type() );
        return resized;
    }


    resizer_type m_resizer;

    wrapped_deriv_type m_dxt;
    wrapped_deriv_type m_dxm;
    wrapped_deriv_type m_dxh;
    wrapped_state_type m_x_tmp;


};

	/************ DUCUMENTATION *******/

} // odeint
} // numeric
} // boost





#endif // BOOST_NUMERIC_ODEINT_STEPPER_RUNGE_KUTTA4_IMPLICIT_HPP_INCLUDED
