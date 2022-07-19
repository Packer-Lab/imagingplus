import os

# from deepinterpolation.generic import JsonSaver, ClassLoader


class Deepinterpolation:
    """module for deepinterpolation denoising on an imaging series"""

    def __init__(self):
        self.generator_param = None
        self.inference_param = None
        self.model_path = None
        self.input_path = None
        self.output_path = None

    @classmethod
    def deepinterpolate(cls, input_file: str, output_file: str, model_path: str, generator_param: dict = {},
                        inferrence_param: dict = {}):
        """
        Alternative constructor for Deepinterpolation module. Used to run and create Deepinterpolation module.

        :param input_file:
        :param output_file:
        :param model_path:
        :param generator_param:
        :param inferrence_param:
        :return:
        """
        _arg_generator = generator_param
        _arg_inferrence = inferrence_param

        generator_param = {}
        inferrence_param = {}

        # We are reusing the data generator for training here.
        generator_param["type"] = "generator"
        generator_param["name"] = "SingleTifGenerator"
        generator_param["pre_post_frame"] = 30
        generator_param["pre_post_omission"] = 0
        generator_param[
            "steps_per_epoch"
        ] = -1
        # No steps necessary for inference as epochs are not relevant.
        # -1 deactivate it.

        generator_param["train_path"] = input_file

        generator_param["batch_size"] = 5
        generator_param["start_frame"] = 0
        generator_param["end_frame"] = -1  # -1 to go until the end.
        generator_param[
            "randomize"
        ] = 0
        # This is important to keep the order
        # and avoid the randomization used during training

        inferrence_param["type"] = "inferrence"
        inferrence_param["name"] = "core_inferrence"

        # Replace this path to where you stored your model
        inferrence_param[
            "model_path"
        ] = model_path

        # Replace this path to where you want to store your output file
        inferrence_param[
            "output_file"
        ] = output_file

        jobdir = os.path.dirname(output_file)

        if len([*_arg_generator]) > 0:
            print('updating generator..')
            for i, val in _arg_generator.items():
                generator_param[i] = val
        if len([*_arg_inferrence]) > 0:
            print('updating inferrence..')
            for i, val in _arg_inferrence.items():
                inferrence_param[i] = val

        try:
            os.mkdir(jobdir)
        except Exception:
            print("folder already exists")

        path_generator = os.path.join(jobdir, "generator.json")
        json_obj = JsonSaver(generator_param)
        json_obj.save_json(path_generator)

        path_infer = os.path.join(jobdir, "inferrence.json")
        json_obj = JsonSaver(inferrence_param)
        json_obj.save_json(path_infer)

        generator_obj = ClassLoader(path_generator)
        data_generator = generator_obj.find_and_build()(path_generator)

        inferrence_obj = ClassLoader(path_infer)
        inferrence_class = inferrence_obj.find_and_build()(path_infer,
                                                           data_generator)

        # Except this to be slow on a laptop without GPU. Inference needs
        # parallelization to be effective.
        inferrence_class.run()

        deepinterpolated = cls()

        deepinterpolated.output_path = output_file
        deepinterpolated.input_path = input_file
        deepinterpolated.model_path = model_path
        deepinterpolated.inference_param = inferrence_param
        deepinterpolated.generator_param = generator_param

        return deepinterpolated
