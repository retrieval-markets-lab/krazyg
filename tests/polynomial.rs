use blstrs::Scalar;
use group::ff::Field;
use krazyg::polynomial::Polynomial;

#[test]
fn test_divide_with_zero_dividend() {
    // Example from https://en.wikipedia.org/wiki/Polynomial_long_division
    let mut divisor = Polynomial::new_zero_poly(4);
    divisor.coeffs[0] = Scalar::from(0) - Scalar::from(3);
    divisor.coeffs[1] = Scalar::from(1);
    let dividend = Polynomial::new_zero_poly(4);
    let (quotient, remainder) = dividend.divide(&divisor, None);
    let expected_quotient= dividend;
    let expected_remainder= Polynomial::new_zero_poly(4);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_divide_with_remainder() {
    // Example from https://en.wikipedia.org/wiki/Polynomial_long_division
    let mut divisor = Polynomial::new_zero_poly(4);
    divisor.coeffs[0] = Scalar::from(0) - Scalar::from(3);
    divisor.coeffs[1] = Scalar::from(1);
    let dividend = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::from(0) - Scalar::from(2), Scalar::one()]);
    let (quotient, remainder) = dividend.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(3), Scalar::one(), Scalar::one(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::from(5), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_divide_with_no_remainder() {
    // Example p(x) = (x+1)(x+2)
    let mut divisor = Polynomial::new_zero_poly(3);
    divisor.coeffs[0] = Scalar::from(1);
    divisor.coeffs[1] = Scalar::from(1);
    let dividend = Polynomial::new(vec![Scalar::from(2), Scalar::from(3), Scalar::one()]);
    let (quotient, remainder) = dividend.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(2), Scalar::one(),Scalar::zero(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_divide_with_no_remainder_2() {
    // Example p(x) = x^2 - 4
    let mut divisor = Polynomial::new_zero_poly(3);
    divisor.coeffs[0] = Scalar::from(0) - Scalar::from(2);
    divisor.coeffs[1] = Scalar::from(1);
    let dividend = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::one()]);
    let (quotient, remainder) = dividend.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(2), Scalar::one(), Scalar::zero(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_divide_with_2x_minus_1_divisor() {
    // Example dividend = x^2 + x, divisor = 2x + 2
    let mut divisor = Polynomial::new_zero_poly(3);
    divisor.coeffs[0] = Scalar::from(2);
    divisor.coeffs[1] = Scalar::from(2);
    let dividend = Polynomial::new(vec![Scalar::from(0), Scalar::one(), Scalar::one()]);
    let (quotient, remainder) = dividend.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::zero(), Scalar::from(2).invert().unwrap(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_find_degree() {
    let dividend_1 = Polynomial::new(vec![Scalar::one(), Scalar::one(), Scalar::one()]);
    let dividend_2 = Polynomial::new(vec![Scalar::one(), Scalar::one(), Scalar::zero()]);
    let dividend_3 = Polynomial::new(vec![Scalar::one(), Scalar::zero(), Scalar::zero()]);
    let expected_1 = Some(2 as usize);
    let expected_2 = Some(1 as usize);
    let expected_3 = Some(0 as usize);
    let output_1 = dividend_1.degree();
    let output_2 = dividend_2.degree();
    let output_3 = dividend_3.degree();
    assert!(output_1 == expected_1);
    assert!(output_2 == expected_2);
    assert!(output_3 == expected_3);
}
#[test]
fn test_find_degree_when_zero_polynomial_given() {
    let zero_poly = Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    assert!(zero_poly.degree() == None);
}