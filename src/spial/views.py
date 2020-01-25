from django.shortcuts import render

def home_page(request):
    context = {
        "title": "Welcome to Spial", 
    }

    return render(request, "home.html", context)
