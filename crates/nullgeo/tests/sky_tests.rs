use nullgeo::EquirectImage;

#[test]
fn equirect_constructor_round_trips_linear_hdr_values() {
    let (w, h) = (4, 2);
    let data = vec![
        [0.0, 0.5, 1.0],
        [12.5, 0.25, 3.0],
        [100.0, 2.0, 0.1],
        [1.0, 1.0, 1.0],
        [0.75, 8.0, 250.0],
        [4.0, 0.0, 16.0],
        [2.5, 32.0, 0.5],
        [1000.0, 1.5, 64.0],
    ];

    let image = EquirectImage::new(w, h, data.clone()).unwrap();
    assert_eq!(image.dimensions(), (w, h));

    for j in 0..h {
        for i in 0..w {
            let expected = data[j * w + i];
            assert_eq!(image.texel(i, j), expected);

            let u = (i as f64 + 0.5) / w as f64;
            let v = (j as f64 + 0.5) / h as f64;
            assert_eq!(image.sample(u, v), expected, "texel center ({i},{j})");
        }
    }
}
