# Copyright 2022 by Owen Bulka.
# All rights reserved.
# This file is released under the "MIT License Agreement".
# Please see the LICENSE.md file that should have been included as part
# of this package.
import logging

import nuke

_LOGGER = logging.getLogger("blink_fog")


class KnobChangedCallbacks(dict):
    def register(self, knob_name):
        def decorated(method):
            self[knob_name] = method
            return method
        return decorated

    def register_multiple(self, knob_names):
        def decorated(method):
            for knob_name in knob_names:
                self[knob_name] = method
            return method
        return decorated


class KnobManager(object):
    _knob_changed_callbacks = KnobChangedCallbacks()

    def __init__(self):
        self._node = nuke.thisNode()
        self._knob = nuke.thisKnob()

    def handle_node_created(self):
        self._input_changed()

    def handle_knob_changed(self):
        knob_name = self._knob.name()
        try:
            type(self)._knob_changed_callbacks.get(knob_name)(self)
        except TypeError:
            _LOGGER.debug("No callbacks for knob: %s", knob_name)

    @_knob_changed_callbacks.register("inputChange")
    def _input_changed(self):
        pass


class BlinkFogKnobManager(KnobManager):
    _knob_changed_callbacks = KnobChangedCallbacks(KnobManager._knob_changed_callbacks)

    _deep_blender = "deep_blender"
    _input_name = "Input1"
    _deep_merge = "DeepMerge"
    _time_warp = "TimeWarp"
    _expression = "Expression"
    _alpha_expression = "ExprAlpha"
    _deep_from_image = "DeepFromImage"
    _holdout_mode = "holdout_mode"
    _samples_per_ray = "samples_per_ray"
    _rays_per_pixel = "rays_per_pixel"
    _output_deep = "output_deep"

    def __init__(self):
        super(BlinkFogKnobManager, self).__init__()

    def _add_deep_samples(self, samples_to_add, current_deep_samples):
        even = current_deep_samples % 2 == 0

        input_pos = nuke.toNode(self._input_name).knob("xpos").value()
        time_warp_y_pos = nuke.toNode(self._time_warp + "0").knob("ypos").value()
        node_spacing_x = 150
        node_spacing_y = 50

        for sample in range(samples_to_add):
            current_sample = current_deep_samples + sample
            x_pos = input_pos + current_sample * node_spacing_x
            if even:
                time_warp = nuke.nodes.TimeWarp()
                time_warp.setInput(0, nuke.toNode(self._input_name))
                time_warp.knob("name").setValue(self._time_warp + "{}".format(current_sample))
                time_warp.knob("xpos").setValue(x_pos)
                time_warp.knob("ypos").setValue(time_warp_y_pos)
                time_warp.knob("filter").setValue("none")
                time_warp.knob("lookup").setExpression("frame + (1/parent.samples_per_ray*{})".format(current_sample))

                expression = nuke.nodes.Expression()
                expression.setInput(0, time_warp)
                expression.knob("name").setValue(self._expression + "{}".format(current_sample))
                expression.knob("xpos").setValue(x_pos)
                expression.knob("ypos").setValue(time_warp_y_pos + node_spacing_y)
                expression.knob("expr0").setValue("r")
                expression.knob("expr1").setValue("r")
                expression.knob("expr2").setValue("r")
                expression.knob("expr3").setValue("1/g")
                expression.knob("channel3").setValue("depth")

                alpha_expression = nuke.nodes.Expression()
                alpha_expression.setInput(0, expression)
                alpha_expression.knob("name").setValue(self._alpha_expression + "{}".format(current_sample))
                alpha_expression.knob("xpos").setValue(x_pos)
                alpha_expression.knob("ypos").setValue(time_warp_y_pos + 2 * node_spacing_y)
                alpha_expression.knob("expr3").setValue("r")

                deep_from_image = nuke.nodes.DeepFromImage()
                deep_from_image.setInput(0, alpha_expression)
                deep_from_image.knob("name").setValue(self._deep_from_image + "{}".format(current_sample))
                deep_from_image.knob("xpos").setValue(x_pos)
                deep_from_image.knob("ypos").setValue(time_warp_y_pos + 3 * node_spacing_y)
                deep_from_image.knob("premult").setValue(True)

                deep_merge = nuke.toNode(self._deep_merge)
                deep_merge.setInput(current_sample, deep_from_image)

            else:
                time_warp = nuke.toNode(self._time_warp + "{}".format(current_sample - 1))

                expression = nuke.nodes.Expression()
                expression.setInput(0, time_warp)
                expression.knob("name").setValue(self._expression + "{}".format(current_sample))
                expression.knob("xpos").setValue(x_pos)
                expression.knob("ypos").setValue(time_warp_y_pos + node_spacing_y)
                expression.knob("expr0").setValue("b")
                expression.knob("expr1").setValue("b")
                expression.knob("expr2").setValue("b")
                expression.knob("expr3").setValue("1/a")
                expression.knob("channel3").setValue("depth")

                alpha_expression = nuke.nodes.Expression()
                alpha_expression.setInput(0, expression)
                alpha_expression.knob("name").setValue(self._alpha_expression + "{}".format(current_sample))
                alpha_expression.knob("xpos").setValue(x_pos)
                alpha_expression.knob("ypos").setValue(time_warp_y_pos + 2 * node_spacing_y)
                alpha_expression.knob("expr3").setValue("r")

                deep_from_image = nuke.nodes.DeepFromImage()
                deep_from_image.setInput(0, alpha_expression)
                deep_from_image.knob("name").setValue(self._deep_from_image + "{}".format(current_sample))
                deep_from_image.knob("xpos").setValue(x_pos)
                deep_from_image.knob("ypos").setValue(time_warp_y_pos + 3 * node_spacing_y)
                deep_from_image.knob("premult").setValue(True)

                deep_merge = nuke.toNode(self._deep_merge)
                deep_merge.setInput(current_sample, deep_from_image)

            even = not even

    def _remove_deep_samples(self, samples_to_remove, current_deep_samples):
        even = current_deep_samples % 2 == 1

        for sample in range(current_deep_samples - 1, current_deep_samples - samples_to_remove - 1, -1):
            if even:
                time_warp = nuke.toNode(self._time_warp + "{}".format(sample))
                time_warp.setInput(0, None)
                nuke.delete(time_warp)

            expression = nuke.toNode(self._expression + "{}".format(sample))
            alpha_expression = nuke.toNode(self._alpha_expression + "{}".format(sample))
            deep_from_image = nuke.toNode(self._deep_from_image + "{}".format(sample))

            nuke.delete(expression)
            nuke.delete(alpha_expression)
            nuke.delete(deep_from_image)

            even = not even

    def update_deep_nodes(self):
        with self._node:
            with nuke.toNode(self._deep_blender):
                current_deep_samples = len(nuke.allNodes(filter="DeepFromImage"))
                samples_per_ray = self._node.knob(self._samples_per_ray).value()
                deep_samples_to_add = int(samples_per_ray - current_deep_samples)

                if deep_samples_to_add > 0:
                    self._add_deep_samples(deep_samples_to_add, current_deep_samples)
                elif deep_samples_to_add < 0 and samples_per_ray > 0:
                    self._remove_deep_samples(abs(deep_samples_to_add), current_deep_samples)

    @_knob_changed_callbacks.register("samples_per_ray")
    def _samples_per_ray_changed(self):
        if self._node.knob(self._holdout_mode).getValue() == 2 or self._node.knob(self._output_deep).value():
            self.update_deep_nodes()

    @_knob_changed_callbacks.register("holdout_mode")
    def _holdout_mode_changed(self):
        with self._node:
            with nuke.toNode(self._deep_blender):
                if self._knob.getValue() < 2 and not self._node.knob(self._output_deep).value():
                    current_deep_samples = len(nuke.allNodes(filter="DeepFromImage"))
                    self._remove_deep_samples(current_deep_samples - 1, current_deep_samples)
                else:
                    self.update_deep_nodes()

    @_knob_changed_callbacks.register("output_deep")
    def _output_deep_changed(self):
        with self._node:
            with nuke.toNode(self._deep_blender):
                if self._knob.value():
                    self.update_deep_nodes()
                elif self._node.knob(self._holdout_mode).getValue() < 2:
                    current_deep_samples = len(nuke.allNodes(filter="DeepFromImage"))
                    self._remove_deep_samples(current_deep_samples - 1, current_deep_samples)
