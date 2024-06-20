---
layout: default
title: "Team"
navtab: "Team"
description: "Team"
group: navigation
navorder: 5
---
{% include JB/setup %}

<style>
    img.photo{
          object-fit: cover;
          border-radius: 50%;
          object-position: 10% 10%; 
          width:150px;
          height:150px;
    }
</style>


<div>
{% include search-form-global.html %}
</div> 

<div class="smalltitle text-left">Research Group Members </div>
<div class="bigspacer"></div>

I'm extremely fortunate to work with several amazing students to whom I serve as primary or co-advisor. <br>

<div class="bigspacer"></div>
<div class="row">
    {% for member in site.categories.team %}
    {% if member.alum == false and member.collaborator == false and member.support == false %}
    <div class="col-sm-3" style="text-align: center">
    {%if member.url%}
    <a href="{{ member.url }}"> <img class="photo" src="{{member.image}}"> </a> <br>
    <div class="head media-heading member-name"><a href="{{ member.url }}" class="off">{{ member.title }}</a></div>  
    <p class="note">{{ member.position }}</p>
    </div>
    {%endif%}
    {%endif%}
    {% endfor %}   
</div>

- [**Jixuan Ni**], Master Student, Computer Science, Rice University.



