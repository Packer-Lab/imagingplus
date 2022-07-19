
def test_run__turboreg(image_registration_fixture):
    ref_path = image_registration_fixture[0]
    data_path = image_registration_fixture[1]

    from packerlabimaging.utils.images import ImportTiff
    from packerlabimaging.processing.turboreg import run__turboreg

    ref_img = ImportTiff(ref_path)[1]
    run__turboreg(stack_path=data_path, ref_img = ref_img)




