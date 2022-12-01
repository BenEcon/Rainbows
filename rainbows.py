from manimlib.imports import *
from math import pi, sqrt
import numpy as np
class Reflection(Scene):
    def construct(self):
        theta = ValueTracker(np.arctan(1/3)) # an increasing number, not, strictly speaking, an angle
        midplane = Line(start = (-6,0,0), end = (6,0,0))
        normal = Line(start = (0,-3,0), end = (0,3,0))
        normal = DashedVMobject(normal)
        inc_start = (-6, 2, 0)
        inc_finish = (0,0, 0)
        in_arrow = Arrow(start = inc_start, end = inc_finish)
        reflec_start = (0,0,0)
        reflec_fin = (-in_arrow.get_start()[0], in_arrow.get_start()[1], 0)
        out_arrow = Arrow(start = reflec_start, end = reflec_fin)
        left_arc = Arc()
        right_arc = Arc()
        right_arc.add_updater(
                lambda m: m.become(
                    Arc(
            radius = 1.5,
            start_angle = PI/2,
            angle = -(PI/2-theta.get_value()), 
            color = BLUE
            )))
        left_arc.add_updater(
            lambda m: m.become(
                Arc(
            radius = 1.5,
            start_angle = PI/2,
            angle = PI/2 - theta.get_value(),
            color = YELLOW
            )))
        label_1 = TexMobject('\\theta_{1}')
        label_2 = TexMobject('\\theta_{2}') 
        label_1.add_updater(
            lambda m: m.move_to(left_arc).shift(.25*UL+1.5*LEFT))
        label_2.add_updater(
            lambda m: m.move_to(right_arc).shift(.25*UR+1.5*RIGHT))
        takeaway = TextMobject("Angle of ", "{incidence }", "= Angle of", "{ reflection}").set_color_by_tex_to_color_map({
            "{incidence }": YELLOW,
            "{ reflection}": BLUE})
        self.add(midplane)
        self.add(normal)
        self.play(GrowArrow(in_arrow))
        in_arrow.add_updater(
            lambda m: 
            m.become(
                # I know what follows is awful but this was the only method I tried that worked.
                Arrow(start = (-(sqrt(40)*np.cos(theta.get_value())), (sqrt(40)*np.sin(theta.get_value())), 0), end = inc_finish)
            ))
        self.play(GrowArrow(out_arrow))
        out_arrow.add_updater(
            lambda m:
            m.become(
                Arrow(start = (reflec_start), end = (-in_arrow.get_start()[0], in_arrow.get_start()[1], 0))
            ))
        self.play(ShowCreation(left_arc), ShowCreation(right_arc))
        self.play(Write(label_1), Write(label_2))
        self.play(theta.increment_value, PI/3- np.arctan(1/3), rate_func = there_and_back, run_time = 2)
        self.play(FadeIn(takeaway.shift(DOWN)))
        self.wait()
        self.remove(left_arc)
        self.wait()
class Snell(Scene):
    def construct(self):
        theta = ValueTracker(np.arctan(1/3)) 
        theta2 = ValueTracker((np.sin(0.75*np.sin(theta.get_value())))) # the refracted angle's value
        midplane = Line(start = (-6,0,0), end = (6,0,0))
        normal = Line(start = (0,-3,0), end = (0,3,0))
        normal = DashedVMobject(normal)
        inc_start = (-6, 2, 0)
        inc_finish = (0,0, 0)
        refrac_start = (0,0,0)
        refrac_fin = (2, -3, 0) # magnitude of this matters but the coordinates really don't...
        in_arrow = Arrow(start = inc_start, end = inc_finish)
        out_arrow = Arrow(start = refrac_start, end = refrac_fin, start_angle = 3*PI/2)
        # ... because the set_angle method effectively changes the coords here, before it is called onscreen
        out_arrow.set_angle(3*PI/2 + np.sin(0.75*np.sin(theta.get_value())))
        left_arc = Arc()
        left_arc.add_updater(
            lambda m: m.become(
                Arc(
            radius = 1.5,
            start_angle = PI/2,
            angle = PI/2 - theta.get_value()
            )))
        right_arc = Arc()
        right_arc.add_updater(
            lambda m: m.become(
                Arc(
            radius = 1.5,
            start_angle = 3*PI/2,
            angle = np.sin(0.75*np.sin(theta.get_value()))
            )))
        label_1 = TexMobject('\\theta_{1}')
        label_2 = TexMobject('\\theta_{2}')
        label_1.add_updater(
            lambda m: m.move_to(left_arc).shift(.165*LEFT+.8*UP))
        label_2.add_updater(
            lambda m: m.move_to(right_arc).shift(.165*RIGHT+.8*DOWN))
        self.add(midplane, normal, in_arrow, left_arc, label_1)
        self.wait()
        self.play(GrowArrow(out_arrow))
        in_arrow.add_updater(
            lambda m: 
            m.become(
                # we use become here because I need to move the arrow's start, which would be fixed if we didn't redefine it like so
                Arrow(start = (-(sqrt(40)*np.cos(theta.get_value())), (sqrt(40)*np.sin(theta.get_value())), 0), end = inc_finish)
                ))
        out_arrow.add_updater(
            # but in this case we can use the set_angle method because all we have to do is rotate it.
            lambda m: m.set_angle(3*PI/2 + np.sin(0.75*np.sin(theta.get_value()))))
        self.play(ShowCreation(right_arc))
        self.play(Write(label_2))
        self.play(theta.increment_value, PI/3- np.arctan(1/3), rate_func = there_and_back, run_time = 2) 
        self.wait()
        self.play(ApplyMethod(in_arrow.rotate, PI))
        self.wait()
class SnellBreakDown(Scene):
    def construct(self):
        initial_snell = TexMobject('n_{1}\\sin\\theta_{1} = n_{2}\\sin\\theta_{2}')
        number_snell = TexMobject('1.00\\sin\\theta_{1} \\approx 1.33\\sin\\theta_{2}')
        frac_snell = TexMobject('\\sin\\theta_{1} = \\frac{4}{3}\\sin\\theta_{2}') 
        self.play(Write(initial_snell))
        self.wait()
        self.play(Transform(initial_snell, number_snell))
        self.wait()
        self.play(Transform(initial_snell, frac_snell))
        self.wait()
class Plot(GraphScene):
    CONFIG = {
        "y_max" : 30,
        "y_min" : 0,
        "x_max" : 180/PI*np.arcsin(3/4),
        "x_min" : 0,
        "y_tick_frequency" : 15, 
        "x_tick_frequency" : 15, 
        "axes_color" : BLUE, 
        "y_labeled_nums": list(np.arange(0, 45, 15)),
        "x_labeled_nums": list(np.arange(0, 60, 15)),
        "x_label_decimal":0,
        "y_label_direction": LEFT,
        "x_label_direction": DOWN,
        "y_label_decimal":0,
        "x_axis_label": "$\\beta$ (degrees)",
        "y_axis_label": "$\\phi$ (degrees)",
    }
    def construct(self):
        dot_x_pos = ValueTracker(180/PI*np.arcsin(np.sqrt(5/12)))
        dot_y_pos = ValueTracker((2*dot_x_pos.get_value()*PI/180 - np.arcsin(4/3*np.sin(dot_x_pos.get_value()*PI/180)))*180/PI)
        self.setup_axes(animate=True)
        graph = self.get_graph(lambda x : (2*x*PI/180 - np.arcsin(4/3*np.sin(x*PI/180)))*180/PI,  
                                    color = GREEN,
                                    )
        self.play(
            ShowCreation(graph),
            run_time = 2
        )
        self.wait()
        vert_line = self.get_vertical_line_to_graph(180/PI*np.arcsin(np.sqrt(5/12)), graph, color=YELLOW).set_opacity(0)
        self.play(GrowArrow(vert_line))
        max_point = Dot(color = RED).move_to(vert_line.get_end())
        self.play(FadeIn(max_point))
        red_arrow = Arrow(start = 3*UP, end = max_point.get_center()+UP*.1, color = RED)
        self.play(FadeIn(red_arrow))
        self.wait()
        self.play(FadeOut(red_arrow))
        vert_line.add_updater(
            lambda m: m.become(
                self.get_vertical_line_to_graph(dot_x_pos.get_value(), graph, color = YELLOW).set_opacity(0)))
        max_point.add_updater(
            lambda m:m.move_to(vert_line.get_end()))
        alpha_param = dot_x_pos.get_value()/self.x_max
        tan_line = TangentLine(vmob = graph, alpha = alpha_param, color = RED, length = 2)
        self.play(GrowArrow(tan_line))
        tan_line.add_updater(
            lambda m: m.become(TangentLine(vmob = graph, alpha = dot_x_pos.get_value()/self.x_max, color = RED, length = 2))
            )
        self.play(dot_x_pos.set_value, 30, rate_func = smooth, run_time = 2)
        self.play(dot_x_pos.set_value, self.x_max-3, rate_func = smooth, run_time = 3)
        self.play(dot_x_pos.set_value, 180/PI*np.arcsin(np.sqrt(5/12)), rate_func = smooth, run_time = 3.25)
        self.wait()
class Eqs(Scene):
    #TODO: add waits EVERYWHERE
    def construct(self):
        dpdb = TexMobject('d\\phi/ d\\beta = 0').shift(2*UP)
        gross = TexMobject('0 = 2 - \\frac{\\frac{4}{3}\\cos\\beta}{\\sqrt{(1- \\frac{16}{9}\\sin^{2}\\beta)}')
        gross2 = TexMobject('2 = \\frac{\\frac{4}{3}\\cos\\beta}{\\sqrt{(1- \\frac{16}{9}\\sin^{2}\\beta)}')
        squared = TexMobject('4 = \\frac{\\frac{16}{9}\\cos^{2}\\beta}{1- \\frac{16}{9}\\sin^{2}\\beta}')
        sub = TexMobject('4 = \\frac{\\frac{16}{9}(1-\\sin^{2}\\beta)}{1- \\frac{16}{9}\\sin^{2}\\beta}}')
        nofrac = TexMobject('4[\\frac{16}{9}(1- \\sin^{2}\\beta)] = (1-\\frac{16}{9}\\sin^{2}\\beta)')
        ellipses = TexMobject('...').shift(1.5*DOWN)
        penult_u = TexMobject('\\sin\\beta_{\\text{max}} = \\sqrt{\\frac{5}{12}}').shift(UP)
        penult_d = TexMobject('\\implies \\beta_{\\text{max}} \\approx 40.2^{\\circ}')
        final_phi = TexMobject('\\implies \\phi_{\\text{max}} \\approx 21.0^{\\circ}').shift(DOWN)
        two_phi = TexMobject('2\\phi_{\\text{max}} \\approx 42^{\\circ}').scale(2).set_color_by_gradient(RED, ORANGE, YELLOW, GREEN, BLUE, PURPLE)
        self.play(Write(dpdb))
        self.wait()
        self.play(Transform(dpdb, gross))
        self.wait()
        self.play(Transform(dpdb, gross2))
        self.wait()
        self.play(Transform(dpdb, squared))
        self.wait()
        self.play(Transform(dpdb, sub))
        self.wait()
        self.play(Transform(dpdb, nofrac))
        self.wait()
        self.play(Write(ellipses), FadeOut(dpdb))
        self.wait()
        self.play(Transform(ellipses, penult_u))
        self.wait()
        self.play(FadeIn(penult_d))
        self.wait()
        self.play(Write(final_phi))
        self.wait()
        self.play(FadeOut(ellipses), FadeOut(penult_d), FadeOut(final_phi))
        # no wait between these two!
        self.play(Write(two_phi))
        self.wait()
# this function and the following class is copied directly from TheoremofBeethoven: 
# https://github.com/Elteoremadebeethoven/MC-TB/blob/master/Tools/Intersection/intersection.py#L13
# His code is open source, but I thought I should mention that. 
# the purpose of putting it here is to make the next few scenes easier
def Range(in_val,end_val,step=1):
    return list(np.arange(in_val,end_val+step,step))
class GetIntersections:
    def get_coord_from_proportion(self,vmob,proportion):
        return vmob.point_from_proportion(proportion)
    def get_points_from_curve(self, vmob, dx=0.005):
        coords = []
        for point in Range(0, 1, dx):
            dot = Dot(self.get_coord_from_proportion(vmob,point))
            coords.append(dot.get_center())
        return coords
    def get_intersections_between_two_vmobs(self, vmob1, vmob2,
                                            tolerance=0.05,
                                            radius_error=0.2,
                                            use_average=True,
                                            use_first_vmob_reference=False):
        coords_1 = self.get_points_from_curve(vmob1)
        coords_2 = self.get_points_from_curve(vmob2)
        intersections = []
        for coord_1 in coords_1:
            for coord_2 in coords_2:
                distance_between_points = get_norm(coord_1 - coord_2)
                if use_average:
                    coord_3 = (coord_2 - coord_1) / 2
                    average_point = coord_1 + coord_3
                else:
                    if use_first_vmob_reference:
                        average_point = coord_1
                    else:
                        average_point = coord_2
                if len(intersections) > 0 and distance_between_points < tolerance:
                    last_intersection=intersections[-1]
                    distance_between_previus_point = get_norm(average_point - last_intersection)
                    if distance_between_previus_point > radius_error:
                        intersections.append(average_point)
                if len(intersections) == 0 and distance_between_points < tolerance:
                    intersections.append(average_point)
        return intersections
class Reprise(Scene, GetIntersections):
    def construct(self):
        # This is where the fun begins
        drop = Circle(radius = 2, color = WHITE).shift(1.5*RIGHT)
        # the y-coords probably won't be true in general
        bisector = DashedLine(start = (drop.get_left()[0], 0,0), end = (drop.get_right()[0], 0, 0))
        in_ray = Line(start = 2.5*UP+2*LEFT, end = (0.5, sqrt(3), 0))
        internal = Line(start = in_ray.get_end(), end = (3.5, 0, 0))
        ray_arrow1 = Arrow(start = in_ray.get_start(), end = in_ray.get_end()).set_width(0.4)
        internal_arrow = Arrow(start = internal.get_start(), end = internal.get_end()).set_width(0.4)
        # a surpise tool that will help us later
        # (more descriptively, the angle from the bisector to the upper normal)
        angle = np.arcsin(in_ray.get_end()[1]/internal.get_length())        
        alpha_param = (PI/2 - angle)/(1*PI)
        normal = TangentLine(vmob = drop, alpha = (alpha_param), color = RED, length = 4).rotate(PI/2)
        normal = DashedLine(start = normal.get_start(), end = drop.get_center()) # i know it's kinda weird but it worked
        horizon = DashedLine(start = (-1+in_ray.get_start()[0], internal.get_start()[1], 0), end = in_ray.get_end())
        # as above, so below
        in_ray_2 = Line(end = 2.5*DOWN+2*LEFT, start = (0.5, -sqrt(3), 0))
        internal_2 = Line(end = in_ray_2.get_start(), start = (3.5, 0, 0))
        ray_arrow2 = Arrow(start = in_ray_2.get_start(), end = in_ray_2.get_end()).set_width(0.4)
        internal_arrow2 = Arrow(start = internal_2.get_start(), end = internal_2.get_end()).set_width(0.4)
        alpha_below = (PI +(2*angle))/(2*PI)
        normal2 = TangentLine(vmob = drop, alpha = (alpha_below), color = RED, length = 4).rotate(PI/2)
        normal2 = DashedLine(start = normal2.get_start(), end = (1.5,0,0)) 
        horizon2 = DashedLine(start = (-1+in_ray_2.get_end()[0], internal_2.get_end()[1], 0), end = in_ray_2.get_start())
        # now to make the angles
        phi_value = np.arctan((2.5-sqrt(3))/2.5)
        beta_value = (normal.get_angle() - internal.get_angle())
        phi_0 = Arc(arc_center = in_ray.get_end(), start_angle = PI, angle = in_ray.get_angle())
        phi_0_label = TexMobject('\\phi').move_to(phi_0).shift(.4*LEFT + .1*UP).scale(0.75)
        phi_0 = VGroup(phi_0, phi_0_label)
        beta_0 = Arc(arc_center = in_ray.get_end(), start_angle = normal.get_angle(), angle = -beta_value)
        beta_0_label = TexMobject('\\beta').move_to(beta_0).scale(0.5).shift(0.3*DR)
        beta_0 = VGroup(beta_0, beta_0_label)
        beta_1 = Arc(arc_center = drop.get_right(), start_angle = PI, angle = -internal.get_angle())
        beta_1_label = TexMobject('\\beta').move_to(beta_1).scale(0.5).shift(0.3*DL+0.1*UP)
        beta_1 = VGroup(beta_1, beta_1_label)
        beta_2 = Arc(arc_center = in_ray_2.get_start(), start_angle = normal2.get_angle(), angle = beta_value)
        beta_2_label = TexMobject('\\beta').move_to(beta_2).scale(0.5).shift(0.3*UR)
        beta_2 = VGroup(beta_2, beta_2_label)
        beta_3 = Arc(arc_center = drop.get_right(), start_angle = PI, angle = internal.get_angle())
        beta_3_label = TexMobject('\\beta').move_to(beta_3).scale(0.5).shift(0.3*UL+0.2*DOWN)
        beta_3 = VGroup(beta_3, beta_3_label)
        phi_1 = Arc(arc_center = in_ray_2.get_start(), start_angle = PI, angle = -in_ray.get_angle())
        phi_1_label = TexMobject('\\phi').move_to(phi_1).shift(.4*LEFT - .1*UP).scale(0.75)
        phi_1 = VGroup(phi_1, phi_1_label)
        red_triangle = Polygon(drop.get_center(), drop.get_right(), in_ray.get_end(), color = RED)
        red_angle_1 = Arc(arc_center = drop.get_center(), start_angle = 0, angle = -2*normal.get_angle(), color = RED, radius = 0.3)
        red_arrow = Arrow(start = 3*UP +3*RIGHT, end = drop.get_center()+.15*UR)
        ra_1_label = TexMobject('180 -2\\beta').move_to(red_arrow.get_start()).shift(1*RIGHT).scale(0.75).set_color(RED)
        # this name is misleading
        red_angle_2 =  Arc(arc_center = drop.get_center(), start_angle = PI, angle = normal.get_angle(), color = BLUE, radius = 0.5)
        ra_2_label = TexMobject('2\\beta').move_to(red_angle_2).shift(0.5*LEFT).scale(0.75).set_color(BLUE)
        big_angle = Arc(arc_center = in_ray.get_end(), start_angle = PI, angle = normal.get_angle(), radius = 1.8, color = BLUE)
        leftmost_label = TexMobject('2\\beta').move_to(big_angle).shift(0.3*UP+.8*LEFT).scale(0.75).set_color(BLUE)
        # if I ever have to make lambda updaters I should probably remove phi_value dependence 
        yellow_angle = Arc(arc_center = in_ray.get_end(), start_angle = PI-phi_value, angle = normal.get_angle()+phi_value, radius = 0.7, color = YELLOW)
        yellow_label = TexMobject('2\\beta - \\phi').move_to(yellow_angle).shift(0.5*UL+0.2*DOWN).scale(0.5).set_color(YELLOW)
        arrow_list = [ray_arrow1, internal_arrow, internal_arrow2, ray_arrow2]
        line_list = [in_ray, internal, internal_2, in_ray_2]
        tupled_list = []
        for i in range(4):
            tupled_list.append((arrow_list[i], line_list[i]))
        snell = TexMobject("\\sin","({2\\beta - \\phi})", "= \\frac{4}{3} \\sin", "({\\beta})").move_to((-4,0,0)).set_color_by_tex_to_color_map({
            "{2\\beta - \\phi}": YELLOW,
            "{\\beta}": RED
            })
        snell2 = TexMobject("({2\\beta - \\phi})", "= \\arcsin\\frac{4}{3}", "( \\sin{\\beta})").move_to((-4,-1.2,0)).set_color_by_tex_to_color_map({
            "{2\\beta - \\phi}": YELLOW,
            "{\\beta}": RED
            })
        snell3 = TexMobject('\\phi = 2\\beta - \\arcsin(\\frac{4}{3} \\sin\\beta)').move_to((-4,0,0))
        self.play(FadeIn(drop))
        for arrow, line in tupled_list:
            self.play(GrowArrow(line), FadeIn(arrow, run_time = 2))
        self.play(GrowArrow(normal))
        self.play(GrowArrow(normal2)) # could maybe combine these two lines
        self.play(FadeIn(bisector))
        self.play(GrowArrow(horizon))
        self.play(GrowArrow(horizon2))
        self.play(FadeIn(phi_0))
        self.play(FadeIn(phi_1))
        up_line = Line(start = in_ray.get_start(), end = in_ray.get_end())
        down_line = Line(start = in_ray_2.get_start(), end = in_ray_2.get_end())
        up_line.scale(7)
        down_line.scale(7)
        intersection = self.get_intersections_between_two_vmobs(up_line, down_line)
        helpful_arrow_up = Line(start = in_ray.get_start(), end = intersection[0], color = GREEN)
        helpful_arrow_down = Line(start = in_ray_2.get_start(), end = intersection[0], color = GREEN)
        helpful_arc = Arc(
            arc_center = intersection[0], 
            start_angle = helpful_arrow_down.get_angle() + PI, 
            angle = helpful_arrow_up.get_angle() - helpful_arrow_down.get_angle(), 
            color = GREEN)
        helpful_label = TexMobject('2\\phi').move_to(helpful_arc.get_center() + 0.5*LEFT)
        self.play(GrowArrow(helpful_arrow_up), GrowArrow(helpful_arrow_down))
        self.play(ShowCreation(helpful_arc))
        self.play(Write(helpful_label))
        self.wait()
class ManyLines(Scene, GetIntersections):
    def construct(self):
        # todo (maybe?): correct the refraction angle 
        drop = Circle(radius = 2, color = WHITE).shift(1.5*RIGHT).rotate(PI)
        self.play(ShowCreation(drop))
        for y in np.linspace(.9, 1.8, num = 6):
            line0 = Line(start = (-2.5, y, 0), end = (0.5*(3-2*sqrt(4 - y**2)),y,0)) # the incident ray
            # these next objects don't show up but are useful for calculation
            invisiline = Line(start = line0.get_end(), end = drop.get_right(), color = BLACK)
            angle = np.arcsin(line0.get_end()[1]/invisiline.get_length()) 
            alpha_param = (PI/2 - angle)/(1*PI)
            invis_normal = TangentLine(vmob = drop, alpha = (alpha_param), color = RED, length = 4).rotate(PI/2)
            new_line = Line(start = line0.get_end(), end = drop.get_right())
            new_line.set_angle(-np.arcsin(3/4 * np.sin(angle)))
            inter_point = self.get_intersections_between_two_vmobs(new_line, drop)
            new_line = Line(start = new_line.get_start(), end = inter_point[1])
            reflec = Line(start = new_line.get_end(), end = drop.get_left()) 
            reflec.set_angle(PI - new_line.get_angle())
            inter_point2 = self.get_intersections_between_two_vmobs(reflec, drop)
            reflec = Line(start = reflec.get_start(), end = inter_point2[-1])
            refrac = Line(start = reflec.get_end())
            refrac.set_angle(reflec.get_angle() + np.arcsin(4/3*np.sin(reflec.get_angle())))
            line_list = [line0, new_line, reflec, refrac]
            for arrow in line_list:
                self.play(GrowArrow(arrow.set_color(RED)))
            self.add(line0)
        self.wait()
class ManyLines2(Scene, GetIntersections):
    def construct(self):
        # this should be exactly the same as the previous class but with different y-values in the for loop
        drop = Circle(radius = 2, color = WHITE).shift(1.5*RIGHT)
        for y in np.linspace(1., 1.4, num = 5):
            line0 = Line(start = (-2.5, y, 0), end = (0.5*(3-2*sqrt(4 - y**2)),y,0))
            invisiline = Line(start = line0.get_end(), end = drop.get_right(), color = BLACK)
            angle = np.arcsin(line0.get_end()[1]/invisiline.get_length())        
            alpha_param = (PI/2 - angle)/(1*PI)
            invis_normal = TangentLine(vmob = drop, alpha = (alpha_param), color = RED, length = 4).rotate(PI/2)
            new_line = Line(start = line0.get_end(), end = drop.get_right())
            new_line.set_angle(-np.arcsin(3/4 * np.sin(angle)))
            inter_point = self.get_intersections_between_two_vmobs(new_line, drop)
            new_line = Line(start = new_line.get_start(), end = inter_point[1])
            reflec = Line(start = new_line.get_end(), end = drop.get_left()) 
            reflec.set_angle(PI - new_line.get_angle())
            inter_point2 = self.get_intersections_between_two_vmobs(reflec, drop)
            reflec = Line(start = reflec.get_start(), end = inter_point2[-1])
            refrac = Line(start = reflec.get_end())
            refrac.set_angle(reflec.get_angle() + np.arcsin(4/3*np.sin(reflec.get_angle())))
            self.add(drop)
            self.play(GrowArrow(line0), GrowArrow(new_line), GrowArrow(reflec), GrowArrow(refrac))
            self.add(line0)
        self.wait()
class PinkFloyd(Scene, GetIntersections):
    def construct(self):
        drop = Circle(radius = 2, color = WHITE).shift(1.5*RIGHT)
        COLOR = {
        'RED':1.33,
        'ORANGE':1.432,
        'YELLOW':1.533,
        'GREEN':1.634,
        'BLUE':1.735,
        'PURPLE':1.837,
        }
        new_line_list = []
        reflec_list = []
        refrac_list = []
        for index in COLOR:
            line0 = Line(start = (-2.5, 1, 0), end = (0.5*(3-2*sqrt(4 - 1**2)),1,0))
            invisiline = Line(start = line0.get_end(), end = drop.get_right(), color = BLACK)
            angle = np.arcsin(line0.get_end()[1]/invisiline.get_length())        
            alpha_param = (PI/2 - angle)/(1*PI)
            invis_normal = TangentLine(vmob = drop, alpha = (alpha_param), color = RED, length = 4).rotate(PI/2)
            new_line = Line(start = line0.get_end(), end = drop.get_right()) #JUST FIX THIS END
            new_line.set_angle(-np.arcsin(1/COLOR.get(index) * np.sin(angle)))
            inter_point = self.get_intersections_between_two_vmobs(new_line, drop)
            new_line = Line(start = new_line.get_start(), end = inter_point[1], color = index)
            reflec = Line(start = new_line.get_end(), end = drop.get_left()) # AND THIS ONE, then everything will work
            reflec.set_angle(PI - new_line.get_angle())
            inter_point2 = self.get_intersections_between_two_vmobs(reflec, drop)
            reflec = Line(start = reflec.get_start(), end = inter_point2[-1], color = index)
            refrac = Line(start = reflec.get_end(), color = index)
            refrac.set_angle(reflec.get_angle() + np.arcsin(1/COLOR.get(index)*np.sin(reflec.get_angle())))
            new_line_list.append(new_line)
            reflec_list.append(reflec)
            refrac_list.append(refrac)
        self.play(FadeIn(drop))
        self.play(GrowArrow(line0))
        self.wait()
        # I tried to find a neater way to do this but it seems impossible since GrowArrow takes only one argument
        self.play(*[GrowArrow(_) for _ in new_line_list])        
        self.wait(1)
        self.play(*[GrowArrow(_) for _ in reflec_list])        
        self.wait(1)
        self.play(*[GrowArrow(_) for _ in refrac_list])        
        self.wait(1)
        self.wait()
class DoubleRainbow(Scene, GetIntersections):
    def construct(self):
        drop = Circle(radius = 2, color = WHITE).shift(1.5*RIGHT)
        self.add(drop)
        inc_down = Line(end = (-2.5, -.75, 0), start = (0.4, -.75, 0))
        inc_down.rotate(-25.4*PI/180).scale(0.7)
        inc_down = Line(start = inc_down.get_end(), end = inc_down.get_start())
        self.wait()
        intersect = self.get_intersections_between_two_vmobs(inc_down, drop)
        inc_down = Line(start = inc_down.get_start(), end = intersect[-1])
        int_down = Line(start = inc_down.get_end(), end = drop.get_right())
        int_down.set_angle(6*PI/180)
        int_down_intersect = self.get_intersections_between_two_vmobs(int_down, drop)
        int_down = Line(start = int_down.get_start(), end = int_down_intersect[-1])
        # as below, so above:
        inc_above = Line(
            start = (inc_down.get_end()[0], -inc_down.get_end()[1], 0),
            end = (inc_down.get_start()[0], -inc_down.get_start()[1], 0)
            )
        int_above = Line(
            end = (int_down.get_start()[0], -int_down.get_start()[1], 0), 
            start = (int_down.get_end()[0], -int_down.get_end()[1], 0)
            )
        int_side = Line(start = int_down.get_end(), end = int_above.get_start())
        arrow_list = [inc_down, int_down, int_side, int_above, inc_above]
        horizon = DashedLine(start = drop.get_left(), end = drop.get_right())
        up_horizon = DashedLine(start = inc_above.get_start(), end = 
            inc_above.get_start()-2*RIGHT)
        down_horizon = DashedLine(start = inc_down.get_end(), end = 
            inc_down.get_end()-2*RIGHT)
        # making/labeling the angles:
        phi_0 = Arc(arc_center = inc_down.get_end(), start_angle = PI, angle = inc_down.get_angle(), radius = 0.6)
        phi_1 = Arc(arc_center = inc_above.get_start(), start_angle = PI, angle = -inc_down.get_angle(), radius = 0.6)
        phi_0_label = TexMobject('\\phi').move_to(phi_0).scale(0.6).shift(LEFT*.25 + 0.1*UP)
        phi_1_label = TexMobject('\\phi').move_to(phi_1).scale(0.6).shift(LEFT*.25 + 0.1*DOWN)
        normal_list = []
        for arrow in arrow_list:
            self.play(GrowArrow(arrow))
            normal = Line(start = arrow.get_start(), end = drop.get_center())
            normal.scale(1.8)
            normal = Line(start = normal.get_start(), end = drop.get_center())
            normal_list.append(DashedLine(start = normal.get_start(), end = normal.get_end()))
        beta_0 = Arc(arc_center = inc_down.get_end(), start_angle = int_down.get_angle(), angle = normal_list[1].get_angle() - int_down.get_angle(), radius = 0.6)
        beta_1 = Arc(arc_center = inc_above.get_start(), start_angle = -int_down.get_angle(), angle = -normal_list[1].get_angle() + int_down.get_angle(), radius = 0.6)
        beta_2 = Arc(arc_center = int_above.get_start(), start_angle = normal_list[3].get_angle(), angle = -(PI/2 + normal_list[3].get_angle()), radius = 0.2)
        beta_3 = Arc(arc_center = int_above.get_start(), start_angle = PI- int_down.get_angle(), angle = 2*PI - np.abs(-normal_list[3].get_angle() + int_above.get_angle()), radius = 0.4)
        perp_bi = Square().move_to((int_side.get_start()[0],0,0)).scale(0.09).shift(.09*UL).rotate(PI/2)
        red_angle = Arc(arc_center = drop.get_center(), angle = PI+(normal_list[3].get_angle()), radius = 0.6, color = RED)
        green_angle = Arc(arc_center = drop.get_center(), start_angle = PI+(normal_list[3].get_angle()), angle = (PI+normal_list[4].get_angle() - (PI+(normal_list[3].get_angle()))), radius = 0.6, color = GREEN) 
        purple_angle = Arc(arc_center = drop.get_center(), start_angle = PI+normal_list[4].get_angle(), angle = -normal_list[4].get_angle(), radius = 0.6, color = PURPLE)
        purple_angle_2 = Arc(arc_center = inc_above.get_start(), start_angle = PI+normal_list[4].get_angle(), angle = -normal_list[4].get_angle(), radius = 0.7, color = PURPLE)
        red_arrow = Arrow(start = drop.get_right()+.5*UP +RIGHT, end = drop.get_center()+.4*RIGHT +.05*UP, color = RED)
        green_arrow = Arrow(start = drop.get_center() + 3*UP +.125*RIGHT, end = drop.get_center()+0.5*UP, color = GREEN)
        red_label = TexMobject('90 - \\beta').move_to(red_arrow.get_start()+0.65*RIGHT).scale(0.6).set_color(RED)
        green_label = TexMobject('180 - 2\\beta').move_to(green_arrow.get_start()+0.3*UP).scale(0.6).set_color(GREEN)
        purple_label = TexMobject('3\\beta - 90').move_to(purple_angle).shift(0.8*LEFT).scale(0.6).set_color(PURPLE)
        purple_label_2 = TexMobject('3\\beta - 90').move_to(purple_angle_2).shift(0.8*LEFT).scale(0.6).set_color(PURPLE)
        beta_0_label = TexMobject('\\beta').move_to(beta_0).scale(0.5).shift(0.25*RIGHT + 0.05*UP)
        beta_1_label = TexMobject('\\beta').move_to(beta_1).scale(0.5).shift(0.25*RIGHT + 0.05*DOWN)
        beta_2_label = TexMobject('\\beta').move_to(beta_2).scale(0.5).shift(0.2*DL)
        beta_3_label = TexMobject('\\beta').move_to(beta_3).scale(0.5).shift(0.2*LEFT+DOWN*.05)
        # ahh, almost done. Now I will create the last equations:
        snell = TexMobject('\\sin(3\\beta - 90^{\\circ} + \\phi) = \\frac{4}{3}\\sin(\\beta)').shift(4*LEFT).scale(0.8)
        phi_eq = TexMobject('\\implies \\phi = 90^{\\circ} - 3\\beta + \\arcsin(\\frac{4}{3}\\beta)').shift(1.5*RIGHT).scale(0.8)
        phi_min = TexMobject('\\phi_{\\text{min}} \\approx 25.4^{\\circ}').shift(DOWN)
        very_last_phi = TexMobject('2\\phi_{\\text{min}} \\approx 51^{\\circ}').shift(DOWN)
        # super dumb but at least it worked:
        very_last_but_bigger = TexMobject('2\\phi_{\\text{min}} \\approx 51^{\\circ}').scale(2)
        self.play(GrowArrow(horizon))
        self.play(GrowArrow(up_horizon), GrowArrow(down_horizon))
        self.play(GrowArrow(normal_list[1]), GrowArrow(normal_list[2]), GrowArrow(normal_list[3]), GrowArrow(normal_list[4]))
        self.play(ShowCreation(phi_0), ShowCreation(beta_0), FadeIn(beta_0_label), FadeIn(phi_0_label))
        self.play(ShowCreation(phi_1), ShowCreation(beta_1), FadeIn(beta_1_label), FadeIn(phi_1_label))
        self.play(ShowCreation(beta_2), FadeIn(beta_2_label))
        self.play(ShowCreation(beta_3), FadeIn(beta_3_label))
        self.play(FadeIn(perp_bi))
        self.play(ShowCreation(red_angle))
        self.play(GrowArrow(red_arrow), Write(red_label))
        self.play(ShowCreation(green_angle))
        self.play(GrowArrow(green_arrow), Write(green_label))
        self.play(ShowCreation(purple_angle), Write(purple_label))
        self.play(ShowCreation(purple_angle_2), Write(purple_label_2))
        self.wait()
        self.play(Write(snell))
        normal_list.pop(0)
        arrow_group = VGroup(*normal_list)
        everything_group = VGroup(
            arrow_group, drop, phi_0, phi_1, horizon, up_horizon, down_horizon, phi_0_label, phi_1_label, beta_0, 
            beta_1, beta_2, beta_3, beta_0_label, beta_1_label, beta_2_label, beta_3_label, perp_bi, inc_down, 
            inc_above, int_down, int_above, red_angle, green_angle, purple_angle, purple_angle_2, red_arrow, 
            green_arrow, red_label, green_label, purple_label, purple_label_2, int_side)
        self.play(FadeOut(everything_group))
        self.play(FadeIn(phi_eq))
        self.play(FadeIn(phi_min))
        self.play(FadeOut(snell), FadeOut(phi_eq))
        self.play(Transform(phi_min, very_last_phi))
        self.play(Transform(phi_min, very_last_but_bigger))
        self.play(ApplyMethod(phi_min.set_color_by_gradient, 
            RED, ORANGE, YELLOW, GREEN, BLUE, PURPLE, 
            RED, ORANGE, YELLOW, GREEN, BLUE, PURPLE))
        self.wait()
class BigPicture(Scene, GetIntersections):
    def construct(self):
        drop = Circle(radius = .8, color = WHITE).shift(0.5*RIGHT +2.8*UP)
        worldline = Line(start = (-7, -2, 0), end = (7,-2,0), color = GREEN)
        flowers = ImageMobject('dorky_flowers.png').scale(0.4).move_to((5, -1.6, 0))
        sticc = ImageMobject('stickman.png').scale(2.3).move_to((-5, -.2, 0))
        COLOR = {
        'RED':1.3,
        'ORANGE':1.35,
        'YELLOW':1.4,
        'GREEN':1.45,
        'BLUE':1.5,
        'PURPLE':1.55,
        }
        new_line_list = []
        reflec_list = []
        refrac_list = []
        for index in COLOR:
            line0 = Line(start = (-3.5, 3.6, 0), end = (8.5, 3.6, 0))
            initial_intersect = self.get_intersections_between_two_vmobs(line0, drop)
            line0 = Line(start = line0.get_start(), end = initial_intersect[0])
            invisiline = Line(start = line0.get_end(), end = drop.get_right(), color = BLACK)
            angle = 40.2*PI/180        
            alpha_param = (PI/2 - angle)/(1*PI)
            invis_normal = TangentLine(vmob = drop, alpha = (alpha_param), color = RED, length = 4).rotate(PI/2)
            new_line = Line(start = line0.get_end(), end = drop.get_right()) 
            new_line.set_angle(-np.arcsin(1/COLOR.get(index) * np.sin(angle)))
            inter_point = self.get_intersections_between_two_vmobs(new_line, drop)
            new_line = Line(start = new_line.get_start(), end = inter_point[1], color = index)
            reflec = Line(start = new_line.get_end(), end = drop.get_left())
            reflec.set_angle(PI - new_line.get_angle())
            inter_point2 = self.get_intersections_between_two_vmobs(reflec, drop)
            reflec = Line(start = reflec.get_start(), end = inter_point2[-1], color = index)
            refrac = Line(start = reflec.get_end(), color = index)
            refrac.set_angle(reflec.get_angle() + np.arcsin(1/COLOR.get(index)*np.sin(reflec.get_angle())))
            refrac.scale(2.4)
            refrac = Line(start = reflec.get_end(), end = refrac.get_end(), color = refrac.get_color())
            new_line_list.append(new_line)
            reflec_list.append(reflec)
            refrac_list.append(refrac)
        self.play(FadeIn(worldline), FadeIn(flowers), FadeIn(sticc))
        self.play(FadeIn(drop))
        self.play(GrowArrow(line0))
        self.wait(1)
        self.play(*[GrowArrow(_) for _ in new_line_list])        
        self.wait(1)
        self.play(*[GrowArrow(_) for _ in reflec_list])        
        self.wait(1)
        self.play(*[GrowArrow(_) for _ in refrac_list])        
        self.wait(1)
        Dashed_1 = DashedLine(start = line0.get_end(), end = (10, 3.6, 0))
        big_invisible_purple = DashedLine(start = refrac_list[2].get_end(), end = (9, 3.6,0))
        big_invisible_purple.set_angle(refrac_list[0].get_angle() + PI)
        disclaimer = TextMobject('*not to scale').move_to((5, 0, 0)).scale(.7)
        self.play(FadeIn(Dashed_1), FadeIn(big_invisible_purple))
        arc_between = Arc(
            arc_center = TOP + RIGHT_SIDE, 
            start_angle = PI + big_invisible_purple.get_angle(), 
            angle = 1.8*big_invisible_purple.get_angle(), 
            radius = 3)
        self.play(ShowCreation(arc_between))
        label_for_arc = TexMobject('42^{\\circ}').move_to(arc_between).shift(LEFT)
        self.play(FadeIn(label_for_arc))
        self.play(FadeIn(disclaimer))
        self.play(FadeOut(disclaimer))
