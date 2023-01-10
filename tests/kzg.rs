use blstrs::Scalar;
use group::ff;
use ff::Field;
use krazyg::{polynomial::Polynomial, kzg::{scalar_powers, KZGPolyCommitmentScheme}};

#[test]
fn test_scalar_powers_1() {
    let expected: Vec<Scalar> = vec![Scalar::one(), Scalar::one(), Scalar::one()];
    let output = scalar_powers(Scalar::one(), 3);
    output.iter().zip(expected.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_scalar_powers_2() {
    let start: Scalar = Scalar::from(2);
    let expected = vec![Scalar::one(), Scalar::from(2), Scalar::from(4)];
    let output = scalar_powers(start, 3);
    output.iter().zip(expected.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_e2e() {
    let kzg_instance = KZGPolyCommitmentScheme::setup(3);
    let poly = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::one()]);
    let commitment = kzg_instance.commit(&poly);
    let witness = kzg_instance.create_witness(&poly, Scalar::from(2));
    let result = kzg_instance.verify(&witness, &commitment);
    assert!(result == true);
}